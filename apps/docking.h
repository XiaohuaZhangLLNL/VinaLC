/* 
 * File:   docking.h
 * Author: zhang30
 *
 * Created on August 14, 2012, 1:53 PM
 */

#ifndef DOCKING_H
#define	DOCKING_H

#ifdef USE_MPI

struct JobInputData{
    bool flexible;
    bool randomize;
    int cpu;
    int exhaustiveness;
    int num_modes;
//    int mc_mult;
    int seed;
    int n[3]; 
    double energy_range;
    double granularity;
    double begin[3];
    double end[3];        
    char ligBuffer[100];
    char ligFile[100000];
    char recBuffer[100];
    char fleBuffer[100];
};

struct JobOutData{       
    char log[1000];
    char poses[100000];
};

int dockjob(JobInputData& jobInput, JobOutData& jobOut);

#endif

struct usage_error : public std::runtime_error {

    usage_error(const std::string & message) : std::runtime_error(message) {
    }
};

#endif	/* DOCKING_H */

