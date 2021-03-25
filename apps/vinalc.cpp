//TODO: remove copyright statement - the license requires referals and License and NOTICE file, not an in-code statement - original authors will still be mentioned
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
//TODO: get rid of all non-needed includes ...
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>
#include <stack>
#include <stdlib.h> // for exit values
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <sys/stat.h> // for kernel retrieved file system info
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
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
#include "coords.h" // add_to_output_container
#include "tokenize.h"


#include <mpi.h>
#include <dock.h>

namespace po = boost::program_options;

//! internal use, only
bool file_exists(const std::string &fname) {
    struct stat finfo;
    int fstat;
    fstat = stat(fname.c_str(), &finfo);
    if (fstat == 0) return true;
    return false;
}

//! simple parsing error to avoid try/catch mangling as in the
//! original VinaLC
void parsingerror(const std::string& msg) {
    std::cerr << std::endl << std::endl
              << "VinaLC was unable to parse your commandline." << std::endl
              << "The parser says: " << std::endl
              << msg << std::endl;
    std::cerr << "Please invoke 'vinalc --help' for detailed information."
              << std::endl;
    exit(EXIT_FAILURE);
}

void inputerror(const std::string& msg) {
    std::cerr << std::endl << std::endl
              << "VinaLC was encountered an Error: " << std::endl
              << msg << std::endl;
}

void helpmessage(po::options_description &option) {
    //TODO: include header to this message
    std::cerr << options << std::endl;
}

void saveStrList(std::string& fileName, std::vector<std::string>& strList){
    std::ifstream inFile;
    try {
        inFile.open(fileName.c_str());
    }
    catch(...){
        std::cout << "Cannot open file" << fileName << std::endl;
    }

    std::string fileLine;
    while(inFile){
        std::getline(inFile, fileLine);
        std::vector<std::string> tokens;
        tokenize(fileLine, tokens);
        if(tokens.size() > 0){
            strList.push_back(tokens[0]);
        }
    }

}

void saveGeoList(std::string& fileName, std::vector<std::vector<double> >& geoList){
    std::ifstream inFile;
    //TODO: replace try/catch with file exist check and abort if not
    try {
        inFile.open(fileName.c_str());
    }
    catch(...){
        std::cout << "Cannot open file" << fileName << std::endl;
    }

    std::string fileLine;
    while(inFile){
        std::getline(inFile, fileLine);
        std::vector<std::string> tokens;
        tokenize(fileLine, tokens);
        if(tokens.size() == 6){
            std::vector<double> geo;
            for(unsigned i=0; i< 6; ++i){
                geo.push_back(atof(tokens[i].c_str()));
            }
            geoList.push_back(geo);
        }
    }
}

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

    #ifdef DEBUG
    std::cout << "Number of tasks= " << nproc
              << " My rank= " << rank << std::endl;
    #endif

    if (rank == 0) {
        #ifdef DEBUG
        std::cout << "Master Node: " << nproc << " My rank= "
                  << rank << std::endl;
        #endif
        std::string recFile;
        std::string fleFile;
        std::string ligFile;
        std::vector<std::string> recList;
        std::vector<std::string> fleList;
        std::vector<std::string> ligList;
        std::vector<std::vector<double> > geoList;
//        std::string recFile;
//        std::string ligFile;
        std::string geoFile;
        bool help;

        options_description inputs("Input");
        inputs.add_options()
                ("recList", po::value<std::string> (&recFile), "receptor list file")
                ("fleList", po::value<std::string> (&fleFile), "flex part receptor list file")
                ("ligList", po::value<std::string> (&ligFile), "ligand list file")
                ("geoList", po::value<std::string> (&geoFile), "receptor geometry file")
                ("exhaustiveness", po::value<int>(&(jobInput.exhaustiveness))->default_value(8), "exhaustiveness (default value 8) of the global search (roughly proportional to time): 1+")
                ("granularity", po::value<double>(&(jobInput.granularity))->default_value(0.375), "the granularity of grids (default value 0.375)")
                ("num_modes", po::value<int>(&jobInput.num_modes)->default_value(9), "maximum number (default value 9) of binding modes to generate")
//                ("mc_mult", value<int>(&jobInput.mc_mult)->default_value(1), "MC step multiplier number (default value 1) [multiply MC steps] ")
                ("seed", po::value<int>(&jobInput.seed), "explicit random seed")
                ("randomize", po::value<bool>(&jobInput.randomize), "Use different random seeds for complex")
                ("energy_range", po::value<fl> (&jobInput.energy_range)->default_value(2.0), "maximum energy difference (default value 2.0) between the best binding mode and the worst one displayed (kcal/mol)")
                ;
        po::options_desctiption runtime("Runtime Parameterization");
        runtime.add_options()
                ("threads,t", po::value<unsigned int>(&jobInput.cpu)->default_value(1), "number of threads per MPI rank (default value: 1)");
        po::options_description info("Information (optional)");
        info.add_options()
                ("help", bool_switch(&help), "display usage summary")
                ;
        po::options_description desc;
        desc.add(inputs).add(info);

        po::variables_map vm;
        try {
            po::store(po::parse_command_line(args, argv, desc), vm);
            po::notify(vm);
        }
        // catching unknown options
        catch (exception_detail::clone_impl<exception_detail::error_info_injector<program_options::unknown_option> > &msg) {
           parsingerror(msg.what());
        }
        // catching missing argument
        catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::invalid_command_line_syntax> > &msg) {
            error(msg.what());
        }
        // invalid option values
        catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::invalid_option_value> > &msg) {
            parsingerror(msg.what());
        }
        // catch multiple occurences of a command line option
        catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::multiple_occurrences> > &msg) {
            string message = " as a command line option";
            parsingerror(msg.what() + message);
        }

        if (vm.count("help") || vm.count("h") {
            helpmessage(desc);
            exit(EXIT_SUCCESS);
        }

        if (vm.count("recList") <= 0) {
            inputerror("Missing receptor list file.");
            helpmessage(desc);
            exit(EXIT_FAILURE);
        } else {
            saveStrList(recFile, recList);
        }

        if (vm.count("fleList") > 0) {
            saveStrList(fleFile, fleList);
        }

        if (vm.count("ligList") <= 0) {
            inputerror("Missing ligand list file.");
            helpmessage(desc);
            exit(EXIT_FAILURE);
        } else {
            saveStrList(ligFile, ligList);
        }

        if (vm.count("geoList") <= 0) {
            inputerror("Missing geometric input data.");
            helpmessage(desc);
            exit(EXIT_FAILURE);
        } else {
            saveGeoList(geoFile, geoList);
        }

        if(geoList.size() != recList.size()){
            inputerror("Receptor and geometry lists are not of equal size.");
            helpmessage(desc);
            exit(EXIT_FAILURE);
        }

        if (vm.count("seed") == 0)
            jobInput.seed = auto_seed();
        //TODO: make this a bool to be parsed
        if (vm.count("randomize") == 0) {
            jobInput.randomize = false;
        } else {
            jobInput.randomize = true;
        }
        if (jobInput.exhaustiveness < 1) {
            inputerror("exhaustiveness must be 1 or greater");
            helpmessage(desc);
            exit(EXIT_FAILURE);
        }
    }
    /* TODO: control whether these catching need to be re-imlemented
    catch (file_error& e) {
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
    */

    // TODO: test whether this can be unsigned int
    int count = 0; // simple counter for job counts
    // TODO: make this string handling stuff more C++alike
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
                        if(count > nproc - 1){
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

                        std::cout << "At Process: " << freeProc << " working on  Ligand: " << ligName.str() << "  receptor: " <<  recList[i] << std::endl;
        //                MPI_Send(recBuffer, 100, MPI_CHAR, freeProc, recTag, MPI_COMM_WORLD);
        //                MPI_Send(geometry, 6, MPI_DOUBLE, freeProc, geoTag, MPI_COMM_WORLD);
                        MPI_Send(&jobInput, sizeof(JobInputData), MPI_CHAR, freeProc, inpTag, MPI_COMM_WORLD);

                    }else {
                        ss << fileLine << std::endl;
                    }
                }
                logFile.flush();
                outFile.flush();
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
        // TODO: output only in verbose mode
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

    #ifdef DEBUG
    time=MPI_Wtime()-time;
    std::cout << "Rank= " << rank << " MPI Wall Time= " << time << std::endl;
    #endif
    MPI_Finalize();
    exit(EXIT_SUCCESS);

}
