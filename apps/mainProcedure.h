/* 
 * File:   mainProcedure.h
 * Author: zhang30
 *
 * Created on September 20, 2012, 4:21 PM
 */

#ifndef MAINPROCEDURE_H
#define	MAINPROCEDURE_H

#include <string>
#include <sstream>

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

//#include "tee.h"
#include "coords.h" // add_to_output_container
#include "tokenize.h"


void doing(int verbosity, const std::string& str, std::stringstream& log);
void done(int verbosity, std::stringstream& log);

path make_path(const std::string& str);

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, 
        const boost::optional<std::string>& flex_name_opt, std::stringstream& ligSS);

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, 
        const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names);

void main_procedure(model& m, const boost::optional<model>& ref, // m is non-const (FIXME?)
        std::stringstream& out_name,
        bool score_only, bool local_only, bool randomize_only, bool no_cache,
        const grid_dims& gd, int exhaustiveness,
        const flv& weights,
        int cpu, int seed, int verbosity, sz num_modes, fl energy_range, std::stringstream& log);

struct usage_error : public std::runtime_error {

    usage_error(const std::string & message) : std::runtime_error(message) {
    }
};

struct options_occurrence {
    bool some;
    bool all;

    options_occurrence() : some(false), all(true) {
    } // convenience

    options_occurrence& operator+=(const options_occurrence & x) {
        some = some || x.some;
        all = all && x.all;
        return *this;
    }
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm, 
        boost::program_options::options_description& d);

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d);

std::string default_output(const std::string& input_name);

#endif	/* MAINPROCEDURE_H */

