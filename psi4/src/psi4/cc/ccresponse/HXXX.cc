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
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {


double HXXX_p1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		      const char *pert_z, int irrep_z, double omega_z) {
    double result = 0.0;
    dpdfile2 X1, mu1, z, z1, l1, mu, lx, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D;
    char lbl[32];
    double Y1_norm;

    global_dpd_->file2_init(&z, PSIF_CC_TMP0, 0, 0, 1, "z_IA");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &z, 0, 0, 1, 0); 
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, 0, 0, 0, "z_IJ");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->contract222(&X1, &z, &z1, 0, 0, 1, 0);
    global_dpd_->file2_close(&z);
    global_dpd_->file2_close(&X1);

    global_dpd_->file2_init(&lx, PSIF_CC_TMP0, 0, 0, 0, "XL_IJ");

    global_dpd_->file2_init(&l1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&l1, &X1, &lx, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&l1);

    result = global_dpd_->file2_dot(&z1, &lx);

    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&lx);


    return result;
}


double HX1X1X1_p1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 Y1, yt, lx;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D, Z2, XW, XL, W;
    char lbl[32];
    double Y1_norm, Y2_norm;    

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "XW(ij,kl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract424(&W, &X1, &XW, 3, 1, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&XW, PSIF_CC_TMP0, rspq, 0, 0, "XW(kl,ij)");
    global_dpd_->buf4_close(&XW);

    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "XL(ij,al)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z(ij,kl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &XL, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&XL);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "XW(kl,ij)");

    result = global_dpd_->buf4_dot(&Z,&XW);

    global_dpd_->buf4_close(&XW);
    global_dpd_->buf4_close(&Z);


    return result; 
}


double HX1X1X1_p2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 Y1, yt, lx;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D, Z2, XW, XL, W;
    char lbl[32];
    double Y1_norm, Y2_norm;

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "XW(ij,kl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract244(&X1, &W, &XW, 1, 3, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&XW);

    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "XL(ij,ka)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&X1, &l2, &XL, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z(ij,kl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &XL, &Z, 1, 3, 0, 1, 0);

    global_dpd_->buf4_close(&XL);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "XW(ij,kl)");

    result = global_dpd_->buf4_dot(&Z,&XW);

    global_dpd_->buf4_close(&XW);
    global_dpd_->buf4_close(&Z);


    //outfile->Printf("\n\tResult G1:  %20.15f\n", result);

    return result;
}



double HX1X1X1_p3(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 Y1, yt, lx;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D, Z2, XW, XL, W;
    char lbl[32];
    double Y1_norm, Y2_norm;

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "XW(aj,bl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract424(&W, &X1, &XW, 3, 1, 1, 1, 0);
//    global_dpd_->contract244(&X1, &W, &XW, 1, 3, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

           Y2_norm = global_dpd_->buf4_dot_self(&XW);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\n\tTODAY XW!!!!! %20.15f\n", Y2_norm);


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 0, 11, 0, 0, "Z(al,kj)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    //global_dpd_->contract244(&X1, &XW, &Z, 1, 2, 0, 1, 0);
    global_dpd_->contract424(&XW, &X1, &Z, 2, 1, 1, 1, 0);
    global_dpd_->buf4_close(&XW);
    global_dpd_->file2_close(&X1);

           Y2_norm = global_dpd_->buf4_dot_self(&Z);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\n\tTODAY Z!!!!! %20.15f\n", Y2_norm);


    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 0, 10, "Z(jk,la)");
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "XL(ij,la)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    //global_dpd_->contract424(&l2, &X1, &XL, 2, 1, 1, 1, 0);
    global_dpd_->contract244(&X1, &l2, &XL, 1, 2, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);


           Y2_norm = global_dpd_->buf4_dot_self(&XL);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\n\tTODAY Z!!!!! %20.15f\n", Y2_norm);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(jk,la)");

    result = global_dpd_->buf4_dot(&XL,&Z);
    //result = global_dpd_->buf4_dot_self(&Z);

    global_dpd_->buf4_close(&XL);
    global_dpd_->buf4_close(&Z);

    outfile->Printf("\n\tResult:  %20.15f\n", result);

    return result;
}


double HX1X1X1_p4(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 Y1, yt, lx;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D, Z2, XW, XL, W;
    char lbl[32];
    double Y1_norm, Y2_norm;

    global_dpd_->buf4_init(&XW, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "XW(aj,lb)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&X1, &W, &XW, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 0, 11, 0, 0, "Z(al,kj)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract424(&XW, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&XW);
    global_dpd_->file2_close(&X1);


    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rspq, 0, 11, "Z(jk,al)");
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_init(&XL, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "XL(ij,al)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract424(&l2, &X1, &XL, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&l2);


    //       Y2_norm = global_dpd_->buf4_dot_self(&XL);
    //       Y2_norm = sqrt(Y2_norm);
    //       outfile->Printf("\n\tTODAY Z!!!!! %20.15f\n", Y2_norm);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(jk,al)");

    result = global_dpd_->buf4_dot(&XL,&Z);

    global_dpd_->buf4_close(&XL);
    global_dpd_->buf4_close(&Z);

    //outfile->Printf("\n\tResult:  %20.15f\n", result);

    return result;
}



double HXXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		     const char *pert_z, int irrep_z, double omega_z) {

    double hyper = 0.0;
    dpdfile2 Y1, yt, lx;
    dpdfile2 X1, mu1, z, z1, l1, mu, lt, xc;
    dpdbuf4 X2, Y2, l2, mu2, z2, Z, D, Z2, W, XW, XL;
    char lbl[32];
    double Y1_norm, Y2_norm;


    /*** <L1(0)|[[[H_bar,X1(A)],X1(B)],X1(C)]|0> ***/ 

    hyper -= HXXX_p1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 
    hyper -= HXXX_p1(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y); 
    hyper -= HXXX_p1(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x); 
    hyper -= HXXX_p1(pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z); 
    hyper -= HXXX_p1(pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);
    hyper -= HXXX_p1(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x); 

    /*** <L2(0)|[[[H_bar,X1(A)],X1(B)],X1(C)]|0> ***/

    hyper += HX1X1X1_p1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += HX1X1X1_p1(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    hyper += HX1X1X1_p1(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);

    hyper += HX1X1X1_p2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper += HX1X1X1_p2(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    hyper += HX1X1X1_p2(pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y, pert_x, irrep_x, omega_x);    

    //hyper -= HX1X1X1_p3(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

    hyper -= HX1X1X1_p4(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);
    hyper -= HX1X1X1_p4(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);
    hyper -= HX1X1X1_p4(pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z, pert_x, irrep_x, omega_x);
    outfile->Printf("\n\tHYPER G1:  %20.15f\n", hyper);


    return hyper;
}

}  // namespace ccresponse
}  // namespace psi
