/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup ccresponse

    Computes CCSD Hyperpolarizability
    Refer to eq. 107 of [Koch:1991:3333] for the general form of quadratic response functions.
 
    Author: Monika Kodrycka      
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <time.h>

namespace psi {
namespace ccresponse {

double YCX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y, 
               const char *pert_z, int irrep_z, double omega_z);
double LCXX(const char *pert_c, int irrep_c, double omega_c, const char *pert_x, int irrep_x, double omega_x,
                     const char *pert_y, int irrep_y, double omega_y); 
double LHXXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                     const char *pert_z, int irrep_z, double omega_z);
double YHXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                     const char *pert_z, int irrep_z, double omega_z);
void quadresp(double *tensor, double A, double B, const char *pert_x, int x_irrep, double omega_x, 
             const char *pert_y, int y_irrep, double omega_y, const char *pert_z, int z_irrep, double omega_z) { 
    double hyper, hyper_YCX;

    /* clear out scratch space */
    for (int j = PSIF_CC_TMP; j <= PSIF_CC_TMP11; j++) {
        psio_close(j, 0);
        psio_open(j, 0);
    }

    hyper = 0.0;

    if ((x_irrep ^ y_irrep ^ z_irrep ) == 0) {
        //if (omega_y != 0.0) { // we assume omega_x = -omega_y 
        timer_on("linear terms");
	//pert A: pert_x
        //pert B: pert_y
        //pert C: pert_z

        //<O|Y1(B)[Abar,X1(C)]|0>
        hyper += YCX(pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x, pert_z, z_irrep, omega_z);

        //<O|Y1(C)[Abar,X1(B)]|0>
        hyper += YCX(pert_z, z_irrep, omega_z, pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);

        //<0|L1(A)[B_bar,X1(C)]|0>
        hyper += YCX(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y, pert_z, z_irrep, omega_z);

        //<0|L1(C)|[B_bar,X1(A)]|0> 
        hyper += YCX(pert_z, z_irrep, omega_z, pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x);

        //<0|L1(A)[C_bar,X1(B)]|0>
        hyper += YCX(pert_x, x_irrep, omega_x, pert_z, z_irrep, omega_z, pert_y, y_irrep, omega_y);

        //<0|L1(B)|[C_bar,X1(A)]|0>
        hyper += YCX(pert_y, y_irrep, omega_y, pert_z, z_irrep, omega_z, pert_x, x_irrep, omega_x);

        timer_off("linear terms");

        hyper += LCXX(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y, pert_z, z_irrep, omega_z);
        hyper += LCXX(pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x, pert_z, z_irrep, omega_z);	
        hyper += LCXX(pert_z, z_irrep, omega_z, pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);       
        hyper += LHXXX(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y, pert_z, z_irrep, omega_z);
        hyper += YHXX(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y, pert_z, z_irrep, omega_z);
        hyper += YHXX(pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x, pert_z, z_irrep, omega_z);
        hyper += YHXX(pert_z, z_irrep, omega_z, pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);

        outfile->Printf("\n\tHyper Final.... %20.15f\n", hyper);	
    }

    *tensor =  hyper; 

}

}  // namespace ccresponse
}  // namespace psi
