/* 
 * File:   mainProcedure.cpp
 * Author: zhang30
 * 
 * Created on September 20, 2012, 4:21 PM
 */

#include "mainProcedure.h"

// which copy
#include <iostream>
#include <fstream>

#include <exception>
#include <stack>
#include <vector> // ligand paths
#include <cmath> // for ceila

using boost::filesystem::path;

path make_path(const std::string& str) {
    //	return path(str, boost::filesystem::native);
    return path(str);
}

void doing(int verbosity, const std::string& str, std::stringstream& log) {
    if (verbosity > 1) {
        log << str << std::string(" ... ");
        log.flush();
    }
}

void done(int verbosity, std::stringstream& log) {
    if (verbosity > 1) {
        log << "done.";
        log << std::endl;
    }
}

std::string default_output(const std::string& input_name) {
    std::string tmp = input_name;
    if (tmp.size() >= 6 && tmp.substr(tmp.size() - 6, 6) == ".pdbqt")
        tmp.resize(tmp.size() - 6); // FIXME?
    return tmp + "_out.pdbqt";
}

void write_all_output(model& m, const output_container& out, sz how_many,
        std::stringstream& output_name,
        const std::vector<std::string>& remarks) {
    if (out.size() < how_many)
        how_many = out.size();
    VINA_CHECK(how_many <= remarks.size());
//    ofile f(make_path(output_name));

    VINA_FOR(i, how_many) {
        m.set(out[i].c);
        m.write_model(output_name, i + 1, remarks[i]); // so that model numbers start with 1
    }
}

void do_randomization(model& m,
        std::stringstream& out_name,
        const vec& corner1, const vec& corner2, int seed, int verbosity, std::stringstream& log) {
    conf init_conf = m.get_initial_conf();
    rng generator(static_cast<rng::result_type> (seed));
    if (verbosity > 1) {
        log << "Using random seed: " << seed;
        log << std::endl;
    }
    const sz attempts = 10000;
    conf best_conf = init_conf;
    fl best_clash_penalty = 0;

    VINA_FOR(i, attempts) {
        conf c = init_conf;
        c.randomize(corner1, corner2, generator);
        m.set(c);
        fl penalty = m.clash_penalty();
        if (i == 0 || penalty < best_clash_penalty) {
            best_conf = c;
            best_clash_penalty = penalty;
        }
    }
    m.set(best_conf);
    if (verbosity > 1) {
        log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
        log << std::endl;
    }
    m.write_structure(out_name);
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc, output_type& out, const vec& cap, sz max_steps = 1000) {
    change g(m.get_size());
    quasi_newton quasi_newton_par;
    quasi_newton_par.max_steps = max_steps;
    const fl slope_orig = nc.slope;

    VINA_FOR(p, 5) {
        nc.slope = 100 * std::pow(10.0, 2.0 * p);
        quasi_newton_par(m, prec, nc, out, g, cap);
        m.set(out.c); // just to be sure
        if (nc.within(m))
            break;
    }
    out.coords = m.get_heavy_atom_movable_coords();
    if (!nc.within(m))
        out.e = max_fl;
    nc.slope = slope_orig;
}

std::string vina_remark(fl e, fl lb, fl ub) {
    std::ostringstream remark;
    remark.setf(std::ios::fixed, std::ios::floatfield);
    remark.setf(std::ios::showpoint);
    remark << "REMARK VINA RESULT: "
            << std::setw(9) << std::setprecision(1) << e
            << "  " << std::setw(9) << std::setprecision(3) << lb
            << "  " << std::setw(9) << std::setprecision(3) << ub
            << '\n';
    return remark.str();
}

output_container remove_redundant(const output_container& in, fl min_rmsd) {
    output_container tmp;
    VINA_FOR_IN(i, in)
    add_to_output_container(tmp, in[i], min_rmsd, in.size());
    return tmp;
}

void do_search(model& m, const boost::optional<model>& ref, const scoring_function& sf, const precalculate& prec, const igrid& ig, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
        std::stringstream& out_name,
        const vec& corner1, const vec& corner2,
        const parallel_mc& par, fl energy_range, sz num_modes,
        int seed, int verbosity, bool score_only, bool local_only, std::stringstream& log, const terms& t, const flv& weights) {
    conf_size s = m.get_size();
    conf c = m.get_initial_conf();
    fl e = max_fl;
    const vec authentic_v(1000, 1000, 1000);
    if (score_only) {
        fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
        naive_non_cache nnc(&prec); // for out of grid issues
        e = m.eval_adjusted(sf, prec, nnc, authentic_v, c, intramolecular_energy);
        log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
        log << std::endl;
        flv term_values = t.evale_robust(m);
        VINA_CHECK(term_values.size() == 5);
        log << "Intermolecular contributions to the terms, before weighting:\n";
        log << std::setprecision(5);
        log << "    gauss 1     : " << term_values[0] << '\n';
        log << "    gauss 2     : " << term_values[1] << '\n';
        log << "    repulsion   : " << term_values[2] << '\n';
        log << "    hydrophobic : " << term_values[3] << '\n';
        log << "    Hydrogen    : " << term_values[4] << '\n';
        VINA_CHECK(weights.size() == term_values.size() + 1);
        fl e2 = 0;
        VINA_FOR_IN(i, term_values)
        e2 += term_values[i] * weights[i];
        e2 = sf.conf_independent(m, e2);
        if (e < 100 && std::abs(e2 - e) > 0.05) {
            log << "WARNING: the individual terms are inconsisent with the\n";
            log << "WARNING: affinity. Consider reporting this as a bug:\n";
            log << "WARNING: http://vina.scripps.edu/manual.html#bugs\n";
        }
    } else if (local_only) {
        output_type out(c, e);
        doing(verbosity, "Performing local search", log);
        refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.evals);
        done(verbosity, log);
        fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out.c);
        e = m.eval_adjusted(sf, prec, nc, authentic_v, out.c, intramolecular_energy);

        log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
        log << std::endl;
        if (!nc.within(m))
            log << "WARNING: not all movable atoms are within the search space\n";

        doing(verbosity, "Writing output", log);
        output_container out_cont;
        out_cont.push_back(new output_type(out));
        std::vector<std::string> remarks(1, vina_remark(e, 0, 0));
        write_all_output(m, out_cont, 1, out_name, remarks); // how_many == 1
        done(verbosity, log);
    } else {
        rng generator(static_cast<rng::result_type> (seed));
        log << "Using random seed: " << seed;
        log << std::endl;
        output_container out_cont;
        doing(verbosity, "Performing search", log);
        par(m, out_cont, prec, ig, prec_widened, ig_widened, corner1, corner2, generator);
        done(verbosity, log);

        doing(verbosity, "Refining results", log);
        VINA_FOR_IN(i, out_cont)
        refine_structure(m, prec, nc, out_cont[i], authentic_v, par.mc.ssd_par.evals);

        if (!out_cont.empty()) {
            out_cont.sort();
            const fl best_mode_intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out_cont[0].c);
            VINA_FOR_IN(i, out_cont)
            if (not_max(out_cont[i].e))
                out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v, out_cont[i].c, best_mode_intramolecular_energy);
            // the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
            out_cont.sort();
        }

        const fl out_min_rmsd = 1;
        out_cont = remove_redundant(out_cont, out_min_rmsd);

        done(verbosity, log);

        log.setf(std::ios::fixed, std::ios::floatfield);
        log.setf(std::ios::showpoint);
        log << '\n';
        log << "mode |   affinity | dist from best mode\n";
        log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
        log << "-----+------------+----------+----------\n";

        model best_mode_model = m;
        if (!out_cont.empty())
            best_mode_model.set(out_cont.front().c);

        sz how_many = 0;
        std::vector<std::string> remarks;

        VINA_FOR_IN(i, out_cont) {
            if (how_many >= num_modes || !not_max(out_cont[i].e) || out_cont[i].e > out_cont[0].e + energy_range) break; // check energy_range sanity FIXME
            ++how_many;
            log << std::setw(4) << i + 1
                    << "    " << std::setw(9) << std::setprecision(1) << out_cont[i].e; // intermolecular_energies[i];
            m.set(out_cont[i].c);
            const model& r = ref ? ref.get() : best_mode_model;
            const fl lb = m.rmsd_lower_bound(r);
            const fl ub = m.rmsd_upper_bound(r);
            log << "  " << std::setw(9) << std::setprecision(3) << lb
                    << "  " << std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures

            remarks.push_back(vina_remark(out_cont[i].e, lb, ub));
            log << std::endl;
        }
        doing(verbosity, "Writing output", log);
        write_all_output(m, out_cont, how_many, out_name, remarks);
        done(verbosity, log);

        if (how_many < 1) {
            log << "WARNING: Could not find any conformations completely within the search space.\n"
                    << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
            log << std::endl;
        }
    }
}

void main_procedure(model& m, const boost::optional<model>& ref, // m is non-const (FIXME?)
        std::stringstream& out_name,
        bool score_only, bool local_only, bool randomize_only, bool no_cache,
        const grid_dims& gd, int exhaustiveness,
        const flv& weights,
        int cpu, int seed, int verbosity, sz num_modes, fl energy_range, std::stringstream& log) {

    doing(verbosity, "Setting up the scoring function", log);

    everything t;
    VINA_CHECK(weights.size() == 6);

    weighted_terms wt(&t, weights);
    precalculate prec(wt);
    const fl left = 0.25;
    const fl right = 0.25;
    precalculate prec_widened(prec);
    prec_widened.widen(left, right);

    done(verbosity, log);

    vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
    vec corner2(gd[0].end, gd[1].end, gd[2].end);

    parallel_mc par;
    sz heuristic = m.num_movable_atoms() + 10 * m.get_size().num_degrees_of_freedom();
    par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
    par.mc.ssd_par.evals = unsigned((25 + m.num_movable_atoms()) / 3);
    par.mc.min_rmsd = 1.0;
    par.mc.num_saved_mins = 20;
    par.mc.hunt_cap = vec(10, 10, 10);
    par.num_tasks = exhaustiveness;
    par.num_threads = cpu;

    const fl slope = 1e6; // FIXME: too large? used to be 100
    if (randomize_only) {
        do_randomization(m, out_name,
                corner1, corner2, seed, verbosity, log);
    } else {
        non_cache nc(m, gd, &prec, slope); // if gd has 0 n's, this will not constrain anything
        non_cache nc_widened(m, gd, &prec_widened, slope); // if gd has 0 n's, this will not constrain anything
        if (no_cache) {
            do_search(m, ref, wt, prec, nc, prec_widened, nc_widened, nc,
                    out_name,
                    corner1, corner2,
                    par, energy_range, num_modes,
                    seed, verbosity, score_only, local_only, log, t, weights);
        } else {
            bool cache_needed = !(score_only || randomize_only || local_only);
            if (cache_needed) doing(verbosity, "Analyzing the binding site", log);
            cache c("scoring_function_version001", gd, slope, atom_type::XS);
            if (cache_needed) c.populate(m, prec, m.get_movable_atom_types(prec.atom_typing_used()));
            if (cache_needed) done(verbosity, log);
            do_search(m, ref, wt, prec, c, prec, c, nc,
                    out_name,
                    corner1, corner2,
                    par, energy_range, num_modes,
                    seed, verbosity, score_only, local_only, log, t, weights);
        }
    }
}


options_occurrence get_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
    options_occurrence tmp;
    VINA_FOR_IN(i, d.options())
    if (vm.count((*d.options()[i]).long_name()))
        tmp.some = true;
    else
        tmp.all = false;
    return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {

    VINA_FOR_IN(i, d.options()) {
        const std::string& str = (*d.options()[i]).long_name();
        if (!vm.count(str))
            std::cerr << "Required parameter --" << str << " is missing!\n";
    }
}

model parse_bundle(const std::string& rigid_name, const boost::optional<std::string>& flex_name_opt, std::stringstream& ligSS) {
    model tmp = (flex_name_opt) ? parse_receptor_pdbqt(make_path(rigid_name), make_path(flex_name_opt.get()))
            : parse_receptor_pdbqt(make_path(rigid_name));
//    VINA_FOR_IN(i, ligand_names)
    tmp.append(parse_ligand_pdbqt(ligSS));
    return tmp;
}

model parse_bundle(const std::string& rigid_name, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
    model tmp = (flex_name_opt) ? parse_receptor_pdbqt(make_path(rigid_name), make_path(flex_name_opt.get()))
            : parse_receptor_pdbqt(make_path(rigid_name));
    VINA_FOR_IN(i, ligand_names)
    tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
    return tmp;
}

model parse_bundle(const std::vector<std::string>& ligand_names) {
    VINA_CHECK(!ligand_names.empty()); // FIXME check elsewhere
    model tmp = parse_ligand_pdbqt(make_path(ligand_names[0]));
    VINA_RANGE(i, 1, ligand_names.size())
    tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
    return tmp;
}

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
    if (rigid_name_opt)
        return parse_bundle(rigid_name_opt.get(), flex_name_opt, ligand_names);
    else
        return parse_bundle(ligand_names);
}

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, const boost::optional<std::string>& flex_name_opt, std::stringstream& ligSS) {
//    if (rigid_name_opt)
        return parse_bundle(rigid_name_opt.get(), flex_name_opt, ligSS);
//    else
//        return parse_bundle(ligand_names);
}

