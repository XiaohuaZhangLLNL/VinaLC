/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

 */
// which copy
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>
#include <stack>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "current_weights.h"
#include "quasi_newton.h"
#include "gzstream.h"
//#include "tee.h"
#include "coords.h" // add_to_output_container
#include "tokenize.h"

#ifdef USE_MPI
#include <mpi.h>
#include <docking.h>
#include <mpiparser.h>
#endif


#ifdef USE_MPI



inline void geometry(JobInputData& jobInput, std::vector<double>& geo){
//    const fl granularity = 0.375;
    vec center(geo[0], geo[1], geo[2]);
    vec span(geo[3], geo[4], geo[5]);

    for(unsigned j=0;j<3; ++j){
        jobInput.n[j]=sz(std::ceil(span[j] / jobInput.granularity));
        fl real_span = jobInput.granularity * jobInput.n[j];
        jobInput.begin[j]=center[j] - real_span / 2;
        jobInput.end[j]=jobInput.begin[j] + real_span;
    }    
}

int main(int argc, char* argv[]) {

    int nproc, rank, rc;

    int jobFlag=1; // 1: doing job,  0: done job
    
    JobInputData jobInput;
    JobOutData jobOut;
    
    MPI_Status status1, status2;
        
    int rankTag=1;
    int jobTag=2;
//    int ligTag=3;
//    int recTag=4;
//    int geoTag=5;
    int inpTag=3;
    int outTag=4;

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        std::cerr << "Error starting MPI program. Terminating.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
//    MPI_Barrier(MPI_COMM_WORLD);
    double time=MPI_Wtime();

    if (nproc < 2) {
        std::cerr << "Error: Total process less than 2" << std::endl;
        return 1;
    }

    std::cout << "Number of tasks= " << nproc << " My rank= " << rank << std::endl;

    if (rank == 0) {
        std::cout << "Master Node: " << nproc << " My rank= " << rank << std::endl;
        std::string recFile;
        std::string fleFile;
        std::string ligFile;      
        std::vector<std::string> recList;
        std::vector<std::string> fleList;
        std::vector<std::string> ligList;
        std::vector<std::vector<double> > geoList;
       
        int success=mpiParser(argc, argv, recFile, fleFile, ligFile, ligList, recList, fleList, geoList, jobInput);
        if(success!=0) {
            std::cerr << "Error: Parser input error" << std::endl;
            return 1;            
        }
        
        unsigned num_cpus = boost::thread::hardware_concurrency();
        if (num_cpus > 0)
            jobInput.cpu = num_cpus;
        else
            jobInput.cpu = 1;    
        
        
        int count=0;
        
        std::string logFName;
        std::string outFName;
        
//        std::ofstream logFile;
//        std::ofstream outFile;
        
        ogzstream logFile;
        ogzstream outFile;
        
        const std::string modStr="MODEL";
        const std::string endStr="ENDMDL";
        
        logFName=recFile+"_"+ligFile+".log.gz";
        logFile.open(logFName.c_str());

        outFName = recFile+"_"+ligFile+".pdbqt.gz";
        outFile.open(outFName.c_str()); 
        
        jobInput.flexible=false;
        if(fleList.size()==recList.size()){
            jobInput.flexible=true;
        }
        
        if(jobInput.randomize){
            srand(unsigned(std::time(NULL)));
        }
        
        for(unsigned i=0; i<recList.size(); ++i){
            std::vector<double> geo=geoList[i];
            geometry(jobInput, geo);
            int ligcount=0; // make the ligand # identical to each receptor.           
            for(unsigned j=0; j<ligList.size(); ++j){
                if(jobInput.randomize){
                    jobInput.seed = rand();
                }
                
                std::ifstream ligFile;
                ligFile.open(ligList[j].c_str());
                
                std::string fileLine;
                std::stringstream ss;
                while(std::getline(ligFile, fileLine)){
                    if(fileLine.compare(0,5, modStr)==0){
                        ss.str(std::string());
                    }else if(fileLine.compare(0,6, endStr)==0){
               
                        ++count;
                        ++ligcount;
                        if(count >nproc-1){
                            MPI_Recv(&jobOut, sizeof(JobOutData), MPI_CHAR, MPI_ANY_SOURCE, outTag, MPI_COMM_WORLD, &status2);
                            logFile << jobOut.log << std::endl;
                            outFile << jobOut.poses << std::endl;
                        }                  
                        int freeProc;
                        MPI_Recv(&freeProc, 1, MPI_INTEGER, MPI_ANY_SOURCE, rankTag, MPI_COMM_WORLD, &status1);
                        MPI_Send(&jobFlag, 1, MPI_INTEGER, freeProc, jobTag, MPI_COMM_WORLD); 
                        // Start to send parameters                        
                        std::stringstream ligName;
                        ligName << "LIGAND " << ligcount;
                        strcpy(jobInput.ligBuffer, ligName.str().c_str());
                        strcpy(jobInput.ligFile, ss.str().c_str());
        //                MPI_Send(ligBuffer, 100, MPI_CHAR, freeProc, ligTag, MPI_COMM_WORLD);  
                        strcpy(jobInput.recBuffer, recList[i].c_str());
                        if(jobInput.flexible){
                                strcpy(jobInput.fleBuffer, fleList[i].c_str());
                        }
        //                MPI_Send(recBuffer, 100, MPI_CHAR, freeProc, recTag, MPI_COMM_WORLD);
        //                MPI_Send(geometry, 6, MPI_DOUBLE, freeProc, geoTag, MPI_COMM_WORLD);
                        MPI_Send(&jobInput, sizeof(JobInputData), MPI_CHAR, freeProc, inpTag, MPI_COMM_WORLD);

                    }else {
                        ss << fileLine << std::endl;
                    }                        
                }  
            }
//            if(i !=recList.size()-1){ // ! Don't close the files if this is last loop for recs.
//                logFile.close();
//                outFile.close();
//            }
        }
        

//        logFile.open(logFName.c_str(), std::ios::app);
//        outFile.open(outFName.c_str(), std::ios::app);
        int nJobs=count;
        int ndata=(nJobs<nproc-1)? nJobs: nproc-1;
        std::cout << "ndata=" << ndata << " nJobs=" << nJobs << std::endl;
    
        for(unsigned i=0; i < ndata; ++i){
            MPI_Recv(&jobOut, sizeof(JobOutData), MPI_CHAR, MPI_ANY_SOURCE, outTag, MPI_COMM_WORLD, &status2);
            logFile << jobOut.log << std::endl;
            outFile << jobOut.poses << std::endl;
        }
        logFile.close();
        outFile.close();
       
        
        for(unsigned i=1; i < nproc; ++i){
            int freeProc;
            MPI_Recv(&freeProc, 1, MPI_INTEGER, MPI_ANY_SOURCE, rankTag, MPI_COMM_WORLD, &status1);
            jobFlag=0;;
            MPI_Send(&jobFlag, 1, MPI_INTEGER, freeProc, jobTag, MPI_COMM_WORLD);            
        }

    } else {
        while (1) {
            MPI_Send(&rank, 1, MPI_INTEGER, 0, rankTag, MPI_COMM_WORLD);
            MPI_Recv(&jobFlag, 20, MPI_CHAR, 0, jobTag, MPI_COMM_WORLD, &status2);
            if (jobFlag==0) {
                break;
            }
            // Receive parameters
//            MPI_Recv(ligBuffer, 100, MPI_CHAR, 0, ligTag, MPI_COMM_WORLD, &status1);
//            MPI_Recv(recBuffer, 100, MPI_CHAR, 0, recTag, MPI_COMM_WORLD, &status1);
//            MPI_Recv(geometry, 6, MPI_DOUBLE, 0, geoTag, MPI_COMM_WORLD, &status1);
            MPI_Recv(&jobInput, sizeof(JobInputData), MPI_CHAR, 0, inpTag, MPI_COMM_WORLD, &status1);
                        
            dockjob(jobInput, jobOut); 
            
            MPI_Send(&jobOut, sizeof(JobOutData), MPI_CHAR, 0, outTag, MPI_COMM_WORLD);
        }
    }

//    MPI_Barrier(MPI_COMM_WORLD);
    time=MPI_Wtime()-time;
    std::cout << "Rank= " << rank << " MPI Wall Time= " << time << std::endl;
    MPI_Finalize();
    return (0);

}


#else

int main(int argc, char* argv[]) {
    using namespace boost::program_options;
    const std::string version_string = "AutoDock Vina 1.1.2 (May 11, 2011)";
    const std::string error_message = "\n\n\
Please contact the author, Dr. Oleg Trott <ot14@columbia.edu>, so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the type of computer system you are running it on,\n\
* all command line options,\n\
* configuration file (if used),\n\
* ligand file as PDBQT,\n\
* receptor file as PDBQT,\n\
* flexible side chains file as PDBQT (if used),\n\
* output file as PDBQT (if any),\n\
* input (if possible),\n\
* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

    const std::string cite_message = "\
#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, Journal of Computational Chemistry 31 (2010)  #\n\
# 455-461                                                       #\n\
#                                                               #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see http://vina.scripps.edu for more information.      #\n\
#################################################################\n";

    try {
        std::string rigid_name, ligand_name, flex_name, config_name, out_name, log_name;
        fl center_x, center_y, center_z, size_x, size_y, size_z;
        int cpu = 0, seed, exhaustiveness, verbosity = 2, num_modes = 9;
        fl energy_range = 2.0;

        // -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
        fl weight_gauss1 = -0.035579;
        fl weight_gauss2 = -0.005156;
        fl weight_repulsion = 0.840245;
        fl weight_hydrophobic = -0.035069;
        fl weight_hydrogen = -0.587439;
        fl weight_rot = 0.05846;
        bool score_only = false, local_only = false, randomize_only = false, help = false, help_advanced = false, version = false; // FIXME

        positional_options_description positional; // remains empty

        options_description inputs("Input");
        inputs.add_options()
                ("receptor", value<std::string > (&rigid_name), "rigid part of the receptor (PDBQT)")
                ("flex", value<std::string > (&flex_name), "flexible side chains, if any (PDBQT)")
                ("ligand", value<std::string > (&ligand_name), "ligand (PDBQT)")
                ;
        //options_description search_area("Search area (required, except with --score_only)");
        options_description search_area("Search space (required)");
        search_area.add_options()
                ("center_x", value<fl > (&center_x), "X coordinate of the center")
                ("center_y", value<fl > (&center_y), "Y coordinate of the center")
                ("center_z", value<fl > (&center_z), "Z coordinate of the center")
                ("size_x", value<fl > (&size_x), "size in the X dimension (Angstroms)")
                ("size_y", value<fl > (&size_y), "size in the Y dimension (Angstroms)")
                ("size_z", value<fl > (&size_z), "size in the Z dimension (Angstroms)")
                ;
        //options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
        options_description outputs("Output (optional)");
        outputs.add_options()
                ("out", value<std::string > (&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
                ("log", value<std::string > (&log_name), "optionally, write log file")
                ;
        options_description advanced("Advanced options (see the manual)");
        advanced.add_options()
                ("score_only", bool_switch(&score_only), "score only - search space can be omitted")
                ("local_only", bool_switch(&local_only), "do local search only")
                ("randomize_only", bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")
                ("weight_gauss1", value<fl > (&weight_gauss1)->default_value(weight_gauss1), "gauss_1 weight")
                ("weight_gauss2", value<fl > (&weight_gauss2)->default_value(weight_gauss2), "gauss_2 weight")
                ("weight_repulsion", value<fl > (&weight_repulsion)->default_value(weight_repulsion), "repulsion weight")
                ("weight_hydrophobic", value<fl > (&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
                ("weight_hydrogen", value<fl > (&weight_hydrogen)->default_value(weight_hydrogen), "Hydrogen bond weight")
                ("weight_rot", value<fl > (&weight_rot)->default_value(weight_rot), "N_rot weight")
                ;
        options_description misc("Misc (optional)");
        misc.add_options()
                ("cpu", value<int>(&cpu), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
                ("seed", value<int>(&seed), "explicit random seed")
                ("exhaustiveness", value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
                ("num_modes", value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
                ("energy_range", value<fl > (&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
                ;
        options_description config("Configuration file (optional)");
        config.add_options()
                ("config", value<std::string > (&config_name), "the above options can be put here")
                ;
        options_description info("Information (optional)");
        info.add_options()
                ("help", bool_switch(&help), "display usage summary")
                ("help_advanced", bool_switch(&help_advanced), "display usage summary with advanced options")
                ("version", bool_switch(&version), "display program version")
                ;
        options_description desc, desc_config, desc_simple;
        desc .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
        desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
        desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

        variables_map vm;
        try {
            //store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
            store(command_line_parser(argc, argv)
                    .options(desc)
                    .style(command_line_style::default_style ^ command_line_style::allow_guessing)
                    .positional(positional)
                    .run(),
                    vm);
            notify(vm);
        } catch (boost::program_options::error& e) {
            std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
            return 1;
        }
        if (vm.count("config")) {
            try {
                path name = make_path(config_name);
                ifile config_stream(name);
                store(parse_config_file(config_stream, desc_config), vm);
                notify(vm);
            } catch (boost::program_options::error& e) {
                std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
                return 1;
            }
        }
        if (help) {
            std::cout << desc_simple << '\n';
            return 0;
        }
        if (help_advanced) {
            std::cout << desc << '\n';
            return 0;
        }
        if (version) {
            std::cout << version_string << '\n';
            return 0;
        }

        bool search_box_needed = !score_only; // randomize_only and local_only still need the search space
        bool output_produced = !score_only;
        bool receptor_needed = !randomize_only;

        if (receptor_needed) {
            if (vm.count("receptor") <= 0) {
                std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
                return 1;
            }
        }
        if (vm.count("ligand") <= 0) {
            std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
            return 1;
        }
        if (cpu < 1)
            cpu = 1;
        if (vm.count("seed") == 0)
            seed = auto_seed();
        if (exhaustiveness < 1)
            throw usage_error("exhaustiveness must be 1 or greater");
        if (num_modes < 1)
            throw usage_error("num_modes must be 1 or greater");
        sz max_modes_sz = static_cast<sz> (num_modes);

        boost::optional<std::string> rigid_name_opt;
        if (vm.count("receptor"))
            rigid_name_opt = rigid_name;

        boost::optional<std::string> flex_name_opt;
        if (vm.count("flex"))
            flex_name_opt = flex_name;

        if (vm.count("flex") && !vm.count("receptor"))
            throw usage_error("Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

        std::stringstream log;
        if (vm.count("log") > 0)
            log.init(log_name);

        if (search_box_needed) {
            options_occurrence oo = get_occurrence(vm, search_area);
            if (!oo.all) {
                check_occurrence(vm, search_area);
                std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
                return 1;
            }
            if (size_x <= 0 || size_y <= 0 || size_z <= 0)
                throw usage_error("Search space dimensions should be positive");
        }

        log << cite_message << '\n';

        if (search_box_needed && size_x * size_y * size_z > 27e3) {
            log << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
        }

        if (output_produced) { // FIXME
            if (!vm.count("out")) {
                out_name = default_output(ligand_name);
                log << "Output will be " << out_name << '\n';
            }
        }

        grid_dims gd; // n's = 0 via default c'tor

        flv weights;
        weights.push_back(weight_gauss1);
        weights.push_back(weight_gauss2);
        weights.push_back(weight_repulsion);
        weights.push_back(weight_hydrophobic);
        weights.push_back(weight_hydrogen);
        weights.push_back(5 * weight_rot / 0.1 - 1); // linearly maps onto a different range, internally. see everything.cpp

        if (search_box_needed) {
            const fl granularity = 0.375;
            vec span(size_x, size_y, size_z);
            vec center(center_x, center_y, center_z);

            VINA_FOR_IN(i, gd) {
                gd[i].n = sz(std::ceil(span[i] / granularity));
                fl real_span = granularity * gd[i].n;
                gd[i].begin = center[i] - real_span / 2;
                gd[i].end = gd[i].begin + real_span;
            }
        }
        if (vm.count("cpu") == 0) {
            unsigned num_cpus = boost::thread::hardware_concurrency();
            if (verbosity > 1) {
                if (num_cpus > 0)
                    log << "Detected " << num_cpus << " CPU" << ((num_cpus > 1) ? "s" : "") << '\n';
                else
                    log << "Could not detect the number of CPUs, using 1\n";
            }
            if (num_cpus > 0)
                cpu = num_cpus;
            else
                cpu = 1;
        }
        if (cpu < 1)
            cpu = 1;
        if (verbosity > 1 && exhaustiveness < cpu)
            log << "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

        doing(verbosity, "Reading input", log);

        model m = parse_bundle(rigid_name_opt, flex_name_opt, std::vector<std::string > (1, ligand_name));

        boost::optional<model> ref;
        done(verbosity, log);

        main_procedure(m, ref,
                out_name,
                score_only, local_only, randomize_only, false, // no_cache == false
                gd, exhaustiveness,
                weights,
                cpu, seed, verbosity, max_modes_sz, energy_range, log);
    } catch (file_error& e) {
        std::cerr << "\n\nError: could not open \"" << e.name.string() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
        return 1;
    } catch (boost::filesystem::filesystem_error& e) {
        std::cerr << "\n\nFile system error: " << e.what() << '\n';
        return 1;
    } catch (usage_error& e) {
        std::cerr << "\n\nUsage error: " << e.what() << ".\n";
        return 1;
    } catch (parse_error& e) {
        std::cerr << "\n\nParse error on line " << e.line << " in file \"" << e.file.string() << "\": " << e.reason << '\n';
        return 1;
    } catch (std::bad_alloc&) {
        std::cerr << "\n\nError: insufficient memory!\n";
        return 1;
    }// Errors that shouldn't happen:
    catch (std::exception& e) {
        std::cerr << "\n\nAn error occurred: " << e.what() << ". " << error_message;
        return 1;
    } catch (internal_error& e) {
        std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << error_message;
        return 1;
    } catch (...) {
        std::cerr << "\n\nAn unknown error occurred. " << error_message;
        return 1;
    }
}

#endif
