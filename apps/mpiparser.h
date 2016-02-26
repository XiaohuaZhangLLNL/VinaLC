/* 
 * File:   mpiparser.h
 * Author: zhang30
 *
 * Created on August 14, 2012, 1:58 PM
 */

#ifndef MPIPARSER_H
#define	MPIPARSER_H
#include <string>

#ifdef USE_MPI
#include "dock.h"

int mpiParser(int argc, char* argv[], 
        std::string& recFile,
        std::string& fleFile,
        std::string& ligFile,
        std::vector<std::string>& ligList,
        std::vector<std::string>& recList,
        std::vector<std::string>& fleList,
        std::vector<std::vector<double> >& geoList,
        JobInputData& jobInput);

#endif
#endif	/* MPIPARSER_H */

