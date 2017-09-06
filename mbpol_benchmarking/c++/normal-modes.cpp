#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstddef>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <string>

#include "constants.h"

#include "io-xyz.h"
#include "io-utils.h"
#include "xyz-water-utils.h"

#include "mbpol.h"

namespace {

const double H_mass = 1.0079;
const double O_mass = 15.9949;

static size_t nw = 0;

std::vector< std::vector<double> > nm;
std::vector<double> freq;
std::vector<std::string> elements;

static x2o::mbpol potential;

//----------------------------------------------------------------------------//

void dump_nm_displacements(const int, double*);

//----------------------------------------------------------------------------//

void hessian_fd(const int natm, double* xyz, double* H)
{
    const double eps = 1.0e-5;

    const size_t ndofs = 3*natm;

    std::fill(H, H + ndofs*ndofs, 0.0);

    for (size_t i = 0; i < ndofs; ++i) {
        const double x0 = xyz[i];

        double gp[ndofs];
        double gm[ndofs];

        xyz[i] = x0 + eps;
	potential(nw, xyz, gp);

        xyz[i] = x0 - eps;
	potential(nw, xyz, gm);

        xyz[i] = x0;

        for (size_t j = 0; j < ndofs; ++j) {
            H[i*ndofs + j] += 0.25*(gp[j] - gm[j])/eps;
            H[j*ndofs + i] += 0.25*(gp[j] - gm[j])/eps;
        }
    }

    //
    // kcal/mol -> J/10 (Na x atomic-mass-unit = 0.001 kg; so that
    //                   lengths are in A and time is in ps)

    for (size_t i = 0; i < ndofs; ++i)
        for (size_t j = 0; j < ndofs; ++j)
            H[j*ndofs + i] *= constants::kcal_J/10.0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: normal-modes optimized-crd.xyz" << std::endl;
        return 0;
    }

    std::vector<double> crd;

    try {
        std::ifstream ifs(argv[1]);

        if (!ifs)
            throw std::runtime_error("could not open the XYZ file");

        std::string comment;
        kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    const int natm = elements.size();
    nw = elements.size()/3;

    const int ndofs = 3*natm;
    double* xyz = &(crd[0]);

    {
        double g[ndofs];
	const double E = potential(nw, xyz, g);

        double gg(0);
        for (size_t n = 0; n < ndofs; ++n)
            gg += g[n]*g[n];

        gg = std::sqrt(gg);

        std::cout << argv[1] << " : E = "
                  << E << " , |g| = "
                  << gg << std::endl;
    }

    double* H = new double[ndofs*ndofs];

    hessian_fd(natm, xyz, H);

    double mass[natm];
    for (size_t i = 0; i < 3*nw; ++i)
	mass[i] = (i%3 == 0 ? O_mass : H_mass);

    for (size_t i = 0; i < natm; ++i) {
	const double mi = mass[i];
        for (size_t j = 0; j < natm; ++j) {
	    const double mj = mass[j];

            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    H[(3*i + k)*ndofs + 3*j + l] /= std::sqrt(mi*mj);
        }
    }

#ifdef VERBOSE
    std::ofstream hessian;
    hessian.open("hessian.dat");

    hessian << std::scientific << std::setprecision(4);
    for(size_t i = 0; i < ndofs; ++i){
	for(size_t j = 0; j < ndofs; ++j)
	    hessian << std::setw(12) << H[i*ndofs + j];

	hessian << std::endl;
    }
    hessian.close();
#endif

    gsl_matrix_view H_gsl = gsl_matrix_view_array(H, ndofs, ndofs);

    gsl_vector* eval = gsl_vector_alloc(ndofs);
    gsl_matrix* evec = gsl_matrix_alloc(ndofs, ndofs);

    gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(ndofs);
    int status = gsl_eigen_symmv(&H_gsl.matrix, eval, evec, ws);
    gsl_eigen_symmv_free(ws);

    std::cout << "GSL::symmv status = " << status << std::endl;

    gsl_eigen_symmv_sort(eval, evec,
                         GSL_EIGEN_SORT_VAL_ASC);

    const double radps_to_cm1 = 1.0e10/(2*M_PI*constants::c0);

    // Print the vibrational frequencies

    std::cout.setf( std::ios::fixed, std::ios::floatfield );
    std::cout << std::setprecision(2);
    for (size_t n = 0; n < ndofs; ++n) {
        double eval_n = gsl_vector_get(eval, n);

	std::cout << std::setw(2) << n;
        if (eval_n > 0.0) {
            const double omega_n = std::sqrt(eval_n); // in rad/ps
            std::cout << "  omega =   "
                      << std::setw(15) << omega_n*radps_to_cm1;
        } else {
            const double omega_n = std::sqrt(-eval_n); // in rad/ps
            std::cout << "  omega = I*"
                      << std::setw(15) << omega_n*radps_to_cm1;
        }

        std::cout << " cm-1" << std::endl;
    }

#ifdef VERBOSE
    // Calculate the cartesian displacements and save them in 'nm'

    for (size_t n = 0; n < ndofs; ++n) {
        double eval_n = gsl_vector_get(eval, n);
        gsl_vector_view evec_n = gsl_matrix_column(evec, n);

	const double this_freq = std::sqrt(std::fabs(eval_n))*radps_to_cm1;

	std::vector<double> this_nm;

	if(this_freq < 1.0){ // assume it is just noise
    	    freq.push_back(0.0);
    	    for (size_t i = 0; i < ndofs; ++i)
		this_nm.push_back(0.0);
	}else{
    	    freq.push_back(this_freq);
    	    for (size_t i = 0; i < natm; ++i)
		for (size_t j = 0; j < 3; ++j)
		    this_nm.push_back( gsl_vector_get(&evec_n.vector, 3*i + j)
				         / std::sqrt(mass[i]));

	    double norm(0);
	    for (size_t i = 0; i < this_nm.size(); ++i)
	    	norm += this_nm[i]*this_nm[i];

	    for (size_t i = 0; i < this_nm.size(); ++i)
		this_nm[i] /= std::sqrt(norm);


	}

	nm.push_back(this_nm);

    }

    // print the displacement matrix

    std::ofstream displacements;
    displacements.open("displacements.dat");

    displacements << std::scientific << std::setprecision(4);
    for(size_t i = 0; i < nm.size(); ++i){
	displacements << std::setw(12) << freq[i];
	for(size_t j = 0; j < nm[i].size(); ++j)
	    displacements << std::setw(12) << nm[i][j];

	displacements << std::endl;
    }


    dump_nm_displacements(natm, xyz);
#endif

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

#ifdef VERBOSE
////////////////////////////////////////////////////////////////////////////////

namespace {

void dump_nm_displacements(const int natm, double* xyz)
{
    const int ndofs = 3*natm;

    const size_t ndisplace = 16;
    double displace[ndisplace] = {0.2, 0.15, 0.125, 
	                          0.1, 0.075, 0.05, 
				  0.0375, 0.025, -0.025, 
				  -0.0375, -0.05, -0.075, 
				  -0.1, -0.125, -0.15, 
				  -0.2};

    for(size_t i = 0; i < ndofs; ++i){
	std::string filename = "nm_" + std::to_string(i) + ".xyz";

	std::ofstream animate;
	animate.open(filename.c_str());

	for(size_t n = 0; n < ndisplace; ++n){
	    double d_xyz[ndofs];
	    std::copy(xyz, xyz + ndofs, d_xyz);

	    for(size_t j = 0; j < ndofs; ++j){
//		std::cerr << d_xyz[j] << ' ' << nm[i][j] << ' ' << displace[n];
		d_xyz[j] += nm[i][j]*displace[n];
//		std::cerr << ' ' << d_xyz[j] 
//		          << std::endl;
	    }

	    tools::print_xyz(animate, elements, d_xyz, std::to_string(freq[i]));
	    
	}
    }
}

} // namespace
#endif
