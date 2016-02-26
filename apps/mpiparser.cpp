#include <iostream>
#include <fstream>

#include <sstream>
#include <exception>
#include <stack>
#include <vector> // ligand paths
#include <cmath> // for ceila
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
//#include "tee.h"
#include "coords.h" // add_to_output_container
#include "tokenize.h"


#include <mpi.h>
#include "dock.h"
#include "mpiparser.h"
#include "mainProcedure.h"


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


int mpiParser(int argc, char* argv[], 
        std::string& recFile,
        std::string& fleFile,
        std::string& ligFile,
        std::vector<std::string>& ligList,
        std::vector<std::string>& recList,
        std::vector<std::string>& fleList,
        std::vector<std::vector<double> >& geoList,
        JobInputData& jobInput){
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
//        std::string recFile;
//        std::string ligFile;
        std::string geoFile;
        bool help;
        
        positional_options_description positional; // remains empty
        
        options_description inputs("Input");
        inputs.add_options()
                ("recList", value<std::string > (&recFile), "receptor list file")
                ("fleList", value<std::string > (&fleFile), "flex part receptor list file")
                ("ligList", value<std::string > (&ligFile), "ligand list file")
                ("geoList", value<std::string > (&geoFile), "receptor geometry file")
                ("exhaustiveness", value<int>(&(jobInput.exhaustiveness))->default_value(8), "exhaustiveness (default value 8) of the global search (roughly proportional to time): 1+")
                ("granularity", value<double>(&(jobInput.granularity))->default_value(0.375), "the granularity of grids (default value 0.375)")
                ("num_modes", value<int>(&jobInput.num_modes)->default_value(9), "maximum number (default value 9) of binding modes to generate")
//                ("mc_mult", value<int>(&jobInput.mc_mult)->default_value(1), "MC step multiplier number (default value 1) [multiply MC steps] ")
                ("seed", value<int>(&jobInput.seed), "explicit random seed")
                ("randomize", value<bool>(&jobInput.randomize), "Use different random seeds for complex")
                ("energy_range", value<fl> (&jobInput.energy_range)->default_value(2.0), "maximum energy difference (default value 2.0) between the best binding mode and the worst one displayed (kcal/mol)")
                ;   
        options_description info("Information (optional)");
        info.add_options()
                ("help", bool_switch(&help), "display usage summary")
                ;
        options_description desc;
        desc.add(inputs).add(info);        

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
            std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc << '\n';
            return 1;
        } 
        
        if (help) {
            std::cout << desc << '\n';
            return 0;
        }    
        
        if (vm.count("recList") <= 0) {
            std::cerr << "Missing receptor List file.\n" << "\nCorrect usage:\n" << desc << '\n';
            return 1;
        }else{
            saveStrList(recFile, recList);
        }
        
        if (vm.count("fleList") > 0) {
            saveStrList(fleFile, fleList);
        }
        
        if (vm.count("ligList") <= 0) {
            std::cerr << "Missing ligand List file.\n" << "\nCorrect usage:\n" << desc << '\n';
            return 1;
        }else{
            saveStrList(ligFile, ligList);
        }  

        if (vm.count("geoList") <= 0) {
            std::cerr << "Missing ligand List file.\n" << "\nCorrect usage:\n" << desc << '\n';
            return 1;
        }else{
            saveGeoList(geoFile, geoList);
        }  
        
        if(geoList.size() != recList.size()){
            std::cerr << "Receptor and geometry lists are not equal.\n" ;
            return 1;        
        }
        
        if (jobInput.cpu < 1)
            jobInput.cpu = 1;
        if (vm.count("seed") == 0)
            jobInput.seed = auto_seed();
        if(vm.count("randomize")==0){
            jobInput.randomize=false;
        }else{
            jobInput.randomize=true;
        }     
        if (jobInput.exhaustiveness < 1)
            throw usage_error("exhaustiveness must be 1 or greater");
        if (jobInput.num_modes < 1)
            throw usage_error("num_modes must be 1 or greater");        
        
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
    
    return 0;
}

