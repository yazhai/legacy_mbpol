#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "io-xyz.h"
#include "io-utils.h"
#include "xyz-water-utils.h"

#include "mbpol.h"

#include "xmin.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//
// cmdline arguments |
//-------------------+

static double grms_tol = 1.0e-6;
static int    maxit = 1000;

static const char* output_file = 0;
static const char* initial_xyz = 0;

void show_usage();
void parse_cmdline(int, char**);

//----------------------------------------------------------------------------//
// global variables |
//------------------+

static size_t nw = 0;

static x2o::mbpol potential;

//----------------------------------------------------------------------------//

double grad_func(double* crd, double* grd, int* natm)
{
    return potential(nw, crd, grd);
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
        parse_cmdline(argc, argv);
        assert(initial_xyz);

        std::ifstream ifs(initial_xyz);

        if (!ifs)
            throw std::runtime_error("could not open the XYZ file");

        std::string comment;
        kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    if (!kit::is_water(elements, crd, false)) {
        std::cerr << " ** Error ** : not water?" << std::endl;
        return 1;
    }

    nw = elements.size()/3;

    double* xyz = new double[2*9*nw];
    double* grd = xyz + 9*nw;

    std::copy(crd.begin(), crd.end(), xyz);

//    std::cout << "<+> Model: " << potential.name() << std::endl;

    kit::xmin_opt xo;
    kit::xmin_opt_init(&xo);

    xo.maxiter = maxit;
    xo.grms_tol = grms_tol;
    xo.print_level = 1;

    xo.method = 3; // 1 - PRCG, 2 -- LBGFS, 3 -- TNCG

    int natm = 3*nw;
    double energy, grms;

    energy = grad_func(xyz, grd, 0);
    const double energy0 = energy;

    kit::xmin(grad_func, &natm, xyz, grd, &energy, &grms, &xo);

    double grd2(0); // to double-check
    for (int n = 0; n < natm; ++n)
        grd2 += grd[n]*grd[n];

    const double energy1 = grad_func(xyz, grd, 0);
    std::cout << "<+> E = " << energy0 << " -> " << energy1
              << ", |g| = " << std::sqrt(grd2)
              << std::endl;

    if (output_file) {
        std::cout << "<+> saving the XYZ to '"
                  << output_file << "'" << std::endl;
        std::ofstream ofs(output_file);

        std::copy(xyz, xyz + 9*nw, crd.begin());
        kit::io::save_xyz(ofs, "optimized", elements, crd);
    }

    delete[] xyz;

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_GETOPT_H
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE 1
#  endif // _GNU_SOURCE
#  include <unistd.h>
#  include <getopt.h>
#else
#  include <unistd.h>
#endif // HAVE_GETOPT_H

extern "C" {
    extern char *optarg;
    extern int optind, opterr, optopt;
} // extern "C"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void show_usage()
{
    std::cout <<
    "\n"
    "Usage: optimize [options] file.xyz\n"
    "\n"
    "Options:\n"
    "  -h --help\n"
    "  -m --maxsteps\n"
    "  -g --grms-tol\n"
    "  -o --output\n"
    "\n";
}

//----------------------------------------------------------------------------//

void parse_cmdline(int argc, char** argv)
{
    using kit::io::to_int;
    using kit::io::to_uint;
    using kit::io::to_double;

    int c;

    static const char short_options[] = "hm:g:o:x:y:";

#   ifdef HAVE_GETOPT_H
    static const struct option long_options[] = {
        {"help",              0, 0, 'h'},
        {"maxsteps",          1, 0, 'm'},
        {"grms-tol",          1, 0, 'g'},
        {"output",            1, 0, 'o'},
        {"x2b-data",          1, 0, 'x'},
        {"x3b-data",          1, 0, 'y'},
        {0,                   0, 0, 0}
    };
#   endif // HAVE_GETOPT_H

    while (true) {
#       ifdef HAVE_GETOPT_H
        c = getopt_long(argc, argv, short_options, long_options, 0);
#       else
        c = getopt(argc, argv, short_options);
#       endif // HAVE_GETOPT_H

        if (c == -1)
            break;

        switch (c) {
            case 'h':
                show_usage();
                std::exit(EXIT_SUCCESS);
            case 'm':
                maxit = to_int(optarg);
                break;
            case 'g':
                grms_tol = to_double(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            case '?':
                std::exit(EXIT_FAILURE);
                break;
            default:
                assert(false); // should not be reached
        }
    } // while (true)

    if (optind >= argc) {
        show_usage();
        std::exit(EXIT_SUCCESS);
    }

    initial_xyz = argv[optind++];

    if (optind + 1 < argc) {
        std::cerr << " ** Warning ** : unexpected command-line arguments:";
        for (; optind < argc; ++optind)
            std::cerr << " '" << argv[optind] << "'";
        std::cerr << std::endl;
    }
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////
