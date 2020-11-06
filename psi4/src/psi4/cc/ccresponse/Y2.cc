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
    \ingroup CCRESPONSE

   Computes inhomogenous terms appearing in Y2 equations.
   
   Author: Monika Kodrycka  

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void denom2(dpdbuf4 *Y2, double omega);
void local_filter_T2(dpdbuf4 *T2);

void Y2_inhomogenous_build(const char *pert, int irrep, double omega) {
    dpdfile2 Y1, X1, L1, mu1, lx_ia, xl, FX, WX, z, F, GAE, GMI, t1, Y1new;  
    dpdbuf4 Y2, Y2new, W, L2, X2, D, LX, XL, lx, lx_ijkb, Z, X1W, Z1, Z2, B, test, test2;
    char lbl[32];
    int Gej, Gab, Gij, Ge, Gi, Gj, Ga, nrows, length, E, e, II;
    int Gbm, Gfe, bm, m, Gb, Gm, Gf, M, fe, f, ef, ncols;
    double *Y;
    dpdbuf4 S, A, B_s;
    int i, j, a, b, ij, ab, Gc, C, c, cc;
    int rows_per_bucket, nbuckets, row_start, rows_left, nlinks;
    psio_address next;
    double Y2_norm;

    //a factor of 0.5 because amplitudes are symmetric
   
    sprintf(lbl, "Inhomo Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");

    // ****<O|L1(0)|A_bar|phi^ab_ij>****

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
 
    global_dpd_->file2_mat_init(&mu1);
    global_dpd_->file2_mat_rd(&mu1);
    global_dpd_->file2_mat_init(&L1);
    global_dpd_->file2_mat_rd(&L1);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z is totally symmetric 
        Gij = Gab ^ irrep;
        global_dpd_->buf4_mat_irrep_init(&Z, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep; 
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z.params->rowtot[Gij]; ij++) {
               i = Z.params->roworb[Gej][ij][0];
               j = Z.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z.params->coltot[Gij]; ab++) {
                   a = Z.params->colorb[Gab][ab][0];
                   b = Z.params->colorb[Gab][ab][1];
                   Z.matrix[Gij][ij][ab]  = 2*L1.matrix[Gi][i][a] * mu1.matrix[Gj][j][b]; 
                   Z.matrix[Gij][ij][ab] -=  L1.matrix[Gj][j][a] * mu1.matrix[Gi][i][b];
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z, Gij);
    }

    global_dpd_->file2_mat_close(&mu1);
    global_dpd_->file2_close(&mu1);
    global_dpd_->file2_mat_close(&L1);
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_axpy(&Z, &Y2new, 1);

    // ****<O|L2(0)|A_bar|phi^ab_ij>****

    sprintf(lbl, "%sBAR_AE", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->contract244(&mu1, &L2, &Y2new, 0, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&mu1);
 
    sprintf(lbl, "%sBAR_MI", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract244(&mu1, &L2, &Y2new, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&mu1);
    global_dpd_->buf4_close(&L2);

    // ****<O|L1(0)|[Hbar(0), X1]|phi^ab_ij>****

    global_dpd_->file2_init(&lx_ia, PSIF_CC_OEI, 0, 0, 1, "Lx_ia");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &X1, &lx_ia, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");

    global_dpd_->file2_mat_init(&lx_ia);
    global_dpd_->file2_mat_rd(&lx_ia);
    global_dpd_->file2_mat_init(&L1);
    global_dpd_->file2_mat_rd(&L1);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z is totally symmetric 
        Gij = Gab ^ irrep;
        global_dpd_->buf4_mat_irrep_init(&Z, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep;
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z.params->rowtot[Gij]; ij++) {
               i = Z.params->roworb[Gej][ij][0];
               j = Z.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z.params->coltot[Gij]; ab++) {
                   a = Z.params->colorb[Gab][ab][0];
                   b = Z.params->colorb[Gab][ab][1];
                   Z.matrix[Gij][ij][ab]  = lx_ia.matrix[Gi][i][b] * L1.matrix[Gj][j][a]; 
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z, Gij);
    }

    global_dpd_->file2_mat_close(&lx_ia);
    global_dpd_->file2_close(&lx_ia);
    global_dpd_->file2_mat_close(&L1);
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_axpy(&Z, &Y2new, -1);

    global_dpd_->file2_init(&xl, PSIF_CC_OEI, 0, 1, 1, "xl_eb");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    global_dpd_->contract222(&X1, &L1, &xl, 1, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>"); 
    global_dpd_->contract424(&D, &xl, &Y2new, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&xl);

    global_dpd_->file2_init(&xl, PSIF_CC_OEI, 0, 0, 0, "xl_mi");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    global_dpd_->contract222(&X1, &L1, &xl, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract244(&xl, &D, &Y2new, 0, 0, 0, -1.0, 1);
    global_dpd_->file2_close(&xl);
    global_dpd_->file2_close(&X1);

    global_dpd_->file2_init(&lx_ia, PSIF_CC_OEI, 0, 0, 1, "Lx_ia");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &X1, &lx_ia, 0, 0, 2.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, " Z(ij|ab)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");

    global_dpd_->file2_mat_init(&lx_ia);
    global_dpd_->file2_mat_rd(&lx_ia);
    global_dpd_->file2_mat_init(&L1);
    global_dpd_->file2_mat_rd(&L1);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z is totally symmetric 
        Gij = Gab ^ irrep;
        global_dpd_->buf4_mat_irrep_init(&Z, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep;
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z.params->rowtot[Gij]; ij++) {
               i = Z.params->roworb[Gej][ij][0];
               j = Z.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z.params->coltot[Gij]; ab++) {
                   a = Z.params->colorb[Gab][ab][0];
                   b = Z.params->colorb[Gab][ab][1];
                   Z.matrix[Gij][ij][ab]  = lx_ia.matrix[Gi][i][a] * L1.matrix[Gj][j][b]; 
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z, Gij);
    }

    global_dpd_->file2_mat_close(&lx_ia);
    global_dpd_->file2_close(&lx_ia);
    global_dpd_->file2_mat_close(&L1);
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_axpy(&Z, &Y2new, 1);

    //****<O|L2(0)|[Hbar(0), X1]|phi^ab_ij>****
 
    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert, omega);  
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract244(&X1, &L2, &LX, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract244(&F, &LX, &Y2new, 0, 2, 1, -1.0, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&LX); 

    global_dpd_->file2_init(&FX, PSIF_CC_OEI, 0, 0, 0, "FX_MI");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract222(&X1, &F, &FX, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&FX, &L2, &Y2new, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&FX);
    global_dpd_->buf4_close(&L2);

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert, omega); 
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&LX, PSIF_CC_LR, qpsr, 0, 11, lbl);
    global_dpd_->buf4_close(&LX);

    sprintf(lbl, "LX_%s_jiak (%5.3f)", pert, omega); 
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract444(&LX, &W, &Y2new, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&LX);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 11, 10, 11, 10, 0, "LX (fj,ma)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&X1, &W, &lx, 1, 2, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);    

    global_dpd_->buf4_sort(&lx, PSIF_CC_LR, qsrp, 10, 10, "LX (ja,mf)");
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z (ib,ja)");

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX (ja,mf)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->contract444(&L2, &lx, &Z, 0, 0, -1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prsq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 11, 11, 11, 11, 0, "LX (fi,bm)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract424(&W, &X1, &lx, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&lx, PSIF_CC_LR, qrsp, 10, 10, "LX (ib,mf)");
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z (ib,ja)");

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX (ib,mf)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->contract444(&lx, &L2, &Z, 0, 0, -1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prsq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 5, 10, "WAmEf 2(Am,Ef) - (Am,fE) (AE,mf)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_init(&WX, PSIF_CC_OEI, irrep, 1, 1, "Wx");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 5, 10, 5, 10, 0, "WAmEf 2(Am,Ef) - (Am,fE) (AE,mf)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&W, &X1, &WX, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);  

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    global_dpd_->contract244( &WX, &L2, &Y2new, 0, 2, 1, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&WX);

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 11, 10, 11, 10, 0, "WX (fi,ma)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); //Compute out of core
    global_dpd_->contract244(&X1, &W, &lx, 1, 2, 1.0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&lx, PSIF_CC_LR, qsrp, 10, 10, "WX (ia,mf)");
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z (ia,jb)");

    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "WX (ia,mf)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->contract444(&lx, &L2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1);
    global_dpd_->buf4_close(&Z);

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "LX_%s_iakj (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&LX, PSIF_CC_LR, psqr, 10, 0, lbl);
    global_dpd_->buf4_close(&LX);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 10, 0, "WMnIe (Me,nI)");  
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z2(ib,ja)");
    sprintf(lbl, "LX_%s_iakj (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMnIe (Me,nI)");      
    global_dpd_->contract444(&LX, &W, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LX);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, prsq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "LX_%s_jikb (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&LX, PSIF_CC_LR, prqs, 0, 10, lbl);
    global_dpd_->buf4_close(&LX);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, sqrp, 11, 0, "WMnIe (en,IM)"); 
    global_dpd_->buf4_close(&W);
    
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "Z2 (ib,aj)");
    sprintf(lbl, "LX_%s_jikb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);
   
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 0, 11, 0, 0, "WMnIe (en,IM)"); 
    global_dpd_->contract444(&LX, &W, &Z2, 1, 0, 1, 0);
    global_dpd_->buf4_close(&LX);
    global_dpd_->buf4_close(&W);   

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, psrq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z2);
    
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "Z (ij,kl)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl); 
    global_dpd_->contract424(&W, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa (Ij,bA)");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z2 (ji,bA)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract444(&Z, &L2, &Z2, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, qpsr, 0, 5, "Z(ij,ab)");

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);
   
    sprintf(lbl, "LX_%s_ijka (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl);

    sprintf(lbl, "LX_%s_kija (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&LX, PSIF_CC_LR, rpqs, 0, 10, lbl);
    global_dpd_->buf4_close(&LX);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 0, 10, "2WMnIe - WnMIe (MI,nE)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (jb,ia)");

    sprintf(lbl, "LX_%s_kija (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&LX, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, lbl); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe (MI,nE)");	
    global_dpd_->contract444(&W, &LX, &Z2, 1, 1, 1, 0);   
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LX);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, prqs, 0, 5, "Z(ij,ab)");

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 0, 10, "2WMnIe - WnMIe (MI,nE)");    
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_init(&z, PSIF_CC_OEI, irrep, 0, 0, "z_ij");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe (MI,nE)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->contract422(&W, &X1, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&z, &L2, &Y2new, 1, 0, 0, -1, 1); 
    global_dpd_->file2_close(&z);
    global_dpd_->buf4_close(&L2);  

    // **** <O|L2(0)|[Hbar(0), X2]|phi^ab_ij> ****

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, "Z (ij,kl)");    
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&L2, &X2, &Z, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");  
    global_dpd_->contract444(&Z, &D, &Y2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, "Z (ij,kl)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D (ij|ba)"); 
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&D, &X2, &Z, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa (Ij,bA)");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa (Ij,bA)");
    global_dpd_->contract444(&Z, &L2, &Y2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->contract444(&L2, &X2, &Z, 1, 1, 1, 0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z2 (ib|ja)"); 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ia|jb)");
    global_dpd_->contract444(&Z, &D, &Z2, 0, 0, 1, 0);   
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prsq, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->contract444(&L2, &X2, &Z, 0, 1, 1, 0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 10, 10, 10, 10, 0, "Z (ib|ja)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ja|ib)"); 
    global_dpd_->contract444(&Z, &D, &Z2, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, prsq, 0, 5, "Z2 (ij,ab)");  
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_jbia");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");  
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&L2, &X2, &Z, 1, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");  
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_jbia");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)"); 
    global_dpd_->contract444(&Z, &D, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, rpsq, 0, 5, "Z2 (ij,ab)"); 
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "G_%s_oo (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&D, &X2, &GMI, 0, 0, 1.0, 0.0); 
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z(ji,ba)"); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    global_dpd_->contract424(&L2, &GMI, &Z2, 1, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&GMI);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, qpsr, 0, 5, "Z(ij,ab)");  
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "G_%s_oo (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&L2, &X2, &GMI, 0, 0, 1.0, 0.0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GMI, &Y2new, 1, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);

    sprintf(lbl, "G_%s_vv (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&D, &X2, &GAE, 2, 2, -1.0, 0.0); 

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&GAE, &L2, &Y2new, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&L2);

    sprintf(lbl, "G_%s_vv (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");   
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&L2, &X2, &GAE, 2, 2, -1.0, 0.0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GAE, &Y2new, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Xl (ib|me)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&L2, &X2, &Z, 1, 0, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);    

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z(ib|ja)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)"); 
    global_dpd_->contract444(&Z, &D, &Z2, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prsq, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z(ib|ja)");
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Xl (ib|me)"); 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");  
    global_dpd_->contract444(&D, &Z, &Z2, 0, 0, 2, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prqs, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_close(&Y2new);
}

}  // namespace ccresponse
}  // namespace psi
