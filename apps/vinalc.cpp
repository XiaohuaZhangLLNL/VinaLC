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


#include <mpi.h>
#include <dock.h>
#include <mpiparser.h>



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
        std::string recFile;
        std::string fleFile;
        std::string ligFile;      
        std::vector<std::string> recList;
        std::vector<std::string> fleList;
        std::vector<std::string> ligList;
        std::vector<std::vector<double> > geoList; 
        int success=mpiParser(argc, argv, recFile, fleFile, ligFile, ligList, recList, fleList, geoList, jobInput);
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


