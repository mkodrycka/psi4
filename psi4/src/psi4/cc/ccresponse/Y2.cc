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
    dpdfile2 Y1, X1, L1, mu1, lx_ia, xl, FX, WX, z, F, GAE, GMI, t1, Y1new;  //xl -> replace this variable
    dpdbuf4 Y2, Y2new, W, L2, X2, D, XL, lx, lx_ijkb, Z, X1W, Z1, Z2, B, test, test2;
    char lbl[32];
    int Gej, Gab, Gij, Ge, Gi, Gj, Ga, nrows, length, E, e, II;
    int Gbm, Gfe, bm, m, Gb, Gm, Gf, M, fe, f, ef, ncols;
    double *Y;
    dpdbuf4 S, A, B_s;
    int i, j, a, b, ij, ab, Gc, C, c, cc;
    int rows_per_bucket, nbuckets, row_start, rows_left, nlinks;
    psio_address next;
    double **Y_diag, **B_diag;
    double Y2_norm;

    //Set of homogenous terms 

    //a factor of 0.5 because teh amplitudes are symmetric
   
    //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega); 
    sprintf(lbl, "Inhomo Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    //sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    //global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    //global_dpd_->buf4_scm(&Y2new, 0); //Do I need this??

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");


    // ****<O|L2(0)|A_bar|phi^ab_ij>****

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
                   Z.matrix[Gij][ij][ab]  = 2*L1.matrix[Gi][i][a] * mu1.matrix[Gj][j][b]; // mu1.matrix[Gi][i][b]; //L1.matrix[Gj][j][a];
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
    //global_dpd_->buf4_close(&L2);
 
    sprintf(lbl, "%sBAR_MI", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract244(&mu1, &L2, &Y2new, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&mu1);
    global_dpd_->buf4_close(&L2);


    // ****<O|L1(0)|[Hbar(0), X1]|phi^ab_ij>****

// I do not know how to make these contractions
    //tmp   = ndot('me,mieb->ib', self.x1, self.Loovv)
    //r_y2 -= ndot('ib,ja->ijab', tmp, self.l1)

    global_dpd_->file2_init(&lx_ia, PSIF_CC_OEI, 0, 0, 1, "Lx_ia");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &X1, &lx_ia, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, irrep, 0, 5, 0, 5, 0, "test");
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
                   Z.matrix[Gij][ij][ab]  = lx_ia.matrix[Gi][i][b] * L1.matrix[Gj][j][a]; // mu1.matrix[Gi][i][b]; //L1.matrix[Gj][j][a];
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

//-----------------------------

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

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract244(&xl, &D, &Y2new, 0, 0, 0, -1.0, 1);
    global_dpd_->file2_close(&xl);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&L1);

// Here 

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
                   Z.matrix[Gij][ij][ab]  = lx_ia.matrix[Gi][i][a] * L1.matrix[Gj][j][b]; // mu1.matrix[Gi][i][b]; //L1.matrix[Gj][j][a];
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
  
    global_dpd_->file2_init(&FX, PSIF_CC_OEI, 0, 1, 1, "FX_AE"); 
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract222(&X1, &F, &FX, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&FX, &L2, &Y2new, 0, 2, 1, -1.0, 1);
    global_dpd_->file2_close(&FX);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_init(&FX, PSIF_CC_OEI, 0, 0, 0, "FX_MI");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract222(&X1, &F, &FX, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);


//sort!!!
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, qpsr, 0, 5, "2 LIjAb - LIjBa (jI,bA)");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa (jI,bA)");
    global_dpd_->contract244(&FX, &L2, &Y2new, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&FX);
    global_dpd_->buf4_close(&L2);


    //tmp   = ndot('me,ijef->mijf', self.x1, self.l2)
    //r_y2 -= ndot('mijf,fmba->ijab', tmp, self.Hvovv) 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, rqps, 11, 10, "(2 LIjAb - LIjBa) (aj|ib)"); //Can we reuse it?
    global_dpd_->buf4_close(&L2); 

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "Z (ij,kb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 11, 10, 11, 10, 0, "(2 LIjAb - LIjBa) (aj|ib)");
    global_dpd_->contract244(&X1, &L2, &Z, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&L2);
   
    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rqsp, 0, 11, "Z (ij,bk)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z2 (ij,ba)");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "Z (ij,bk)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract444(&Z, &W, &Z2, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);    

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, pqsr, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");    
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);  
    global_dpd_->buf4_close(&Z2);

//END HERE
//---------------------------------------------
/* It is NOT WORKING!!!
//sort!!
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "(2 LIjAb - LIjBa) (ij|ba)"); //Can we reuse it?
    global_dpd_->buf4_close(&L2); 


    global_dpd_->buf4_init(&lx_ijkb, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, "LX_ijkb"); 
    //global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "(2 LIjAb - LIjBa) (ij|ba)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
outfile->Printf("\n I am here3a");
      global_dpd_->contract244(&X1, &L2, &lx_ijkb, 1, 3, 0, 1.0, 0); //???
    //global_dpd_->contract424(&L2, &X1, &lx_ijkb, 2, 1, 1, 1.0, 0); //???
outfile->Printf("\n I am here4a");

           Y2_norm = global_dpd_->buf4_dot_self(&lx_ijkb);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\tNorm of lx_ijkb begin %20.15f\n", Y2_norm);

//---
    global_dpd_->buf4_sort(&lx_ijkb, PSIF_CC_LR, pqsr, 0, 11, "LX_ijbk");
    global_dpd_->buf4_close(&lx_ijkb);

    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_init(&test, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "test");
    //global_dpd_->buf4_init(&lx_ijkb, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "LX_ijkb");
    global_dpd_->buf4_init(&lx_ijkb, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "LX_ijbk");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (am,fe)"); //Compute this part out of core /I am reusing!
    //gilobal_dpd_->contract444(&lx_ijkb, &W, &Y2new, 0, 1, -1, 1);
    global_dpd_->contract444(&lx_ijkb, &W, &test, 0, 1, -1, 0);

    global_dpd_->buf4_close(&lx_ijkb);
    global_dpd_->buf4_close(&W); 


    Y2_norm = global_dpd_->buf4_dot_self(&test);
    Y2_norm = sqrt(Y2_norm);
    outfile->Printf("\tNorm of test!!!!! begin %20.15f\n\n", Y2_norm);

   global_dpd_->buf4_close(&test);
*/


// It is not working!!!!

/*
//sort!!
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, rpqs, 11, 10, "(2 LIjAb - LIjBa) (ai|jb)");
    global_dpd_->buf4_close(&L2); 

    global_dpd_->buf4_init(&lx_ijkb, PSIF_CC_LR, irrep, 0, 10, 0, 10, 0, "LX_ijkb"); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 11, 10, 11, 10, 0, "(2 LIjAb - LIjBa) (ai|jb)"); 
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
outfile->Printf("\n I am here3a");
    //global_dpd_->contract244(&X1, &L2, &lx_ijkb, 1, 0, 1, 1.0, 0); //???
    global_dpd_->contract424(&L2, &X1, &lx_ijkb, 0, 1, 1, 1.0, 0); //???
outfile->Printf("\n I am here4a");

           Y2_norm = global_dpd_->buf4_dot_self(&lx_ijkb);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\tNorm of lx_ijkb begin %20.15f\n", Y2_norm);

    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);


    //sort lx_ijkb //can we reuse it?
    global_dpd_->buf4_sort(&lx_ijkb, PSIF_CC_LR, qrsp, 0, 11, "LX_ijkb (jk,bi)"); //Can we reuse it?
    global_dpd_->buf4_close(&lx_ijkb);


outfile->Printf("\n I am here5");

    global_dpd_->buf4_init(&test, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "test");
    global_dpd_->buf4_init(&lx_ijkb, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "LX_ijkb (jk,bi)");
     global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (am,fe)"); //Compute this part out of core /I am reusing!
    //gilobal_dpd_->contract444(&lx_ijkb, &W, &Y2new, 0, 1, -1, 1);
    global_dpd_->contract444(&lx_ijkb, &W, &test, 0, 1, -1, 0);

    Y2_norm = global_dpd_->buf4_dot_self(&W);
    Y2_norm = sqrt(Y2_norm);
    outfile->Printf("\tNorm of test W!!!!! begin %20.15f\n", Y2_norm);


    global_dpd_->buf4_close(&lx_ijkb);
    global_dpd_->buf4_close(&W); 


    Y2_norm = global_dpd_->buf4_dot_self(&test);
    Y2_norm = sqrt(Y2_norm);
    outfile->Printf("\tNorm of test!!!!! begin %20.15f\n", Y2_norm);

   global_dpd_->buf4_close(&test);
*/

//------------------STOP------------------------------------------------

//------------------------------
    //tmp   = ndot('me,imbf->eibf', self.x1, self.l2)
    //r_y2 -= ndot('eibf,fjea->ijab', tmp, self.Hvovv)
//sort!!
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, qprs, 0, 5, "2 LIjAb - LIjBa (ji|ab)");
    global_dpd_->buf4_close(&L2);


    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 11, 5, 11, 5, 0, "LX_eibf");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa (ji|ab)");
    global_dpd_->contract244(&X1, &L2, &lx, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&L2);

//sort 
    global_dpd_->buf4_sort(&lx, PSIF_CC_LR, qrsp, 10, 5, "LX_eibf (ib,fe)");
    global_dpd_->buf4_close(&lx);

//sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 5, 10, "WAmEf (AE,mf)");
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z (ib,ja))");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 5, 10, 5, 10, 0, "WAmEf (AE,mf)");
    global_dpd_->buf4_init(&lx, PSIF_CC_LR, irrep, 10, 5, 10, 5, 0, "LX_eibf (ib,fe)");
    global_dpd_->contract444(&lx, &W, &Z, 0, 1, -1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&lx);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prsq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1);
    global_dpd_->buf4_close(&Z);


    //tmp   = ndot('me,jmfa->ejfa', self.x1, self.l2)
    //r_y2 -= ndot('fibe,ejfa->ijab', self.Hvovv, tmp)
    global_dpd_->buf4_init(&X1W, PSIF_CC_HBAR, 0, 11, 11, 11, 11, 0, "X1W");	
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract424(&W, &X1, &X1W, 3, 1, 0.0, 1.0, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_sort(&X1W, PSIF_CC_HBAR, rqps, 11, 11, "Z (bi,fm)");
    global_dpd_->buf4_close(&X1W);   

    //sort  L
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, psrq, 10, 11, "(2 LIjAb - LIjBa) (ib|aj)"); 
    global_dpd_->buf4_close(&L2);
   
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, irrep, 11, 10, 11, 10, 0, "Z (bi,ja)");
    global_dpd_->buf4_init(&X1W, PSIF_CC_HBAR, 0, 11, 11, 11, 11, 0, "Z (bi,fm)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 11, 10, 11, 0, "(2 LIjAb - LIjBa) (ib|aj)");
    global_dpd_->contract444(&X1W, &L2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X1W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, qrsp, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, -1);
    global_dpd_->buf4_close(&Z);


    //tmp   = ndot('me,fmae->fa', self.x1, self.Hvovv, prefactor=2.0)
    //tmp  -= ndot('me,fmea->fa', self.x1, self.Hvovv)
    //r_y2 += ndot('ijfb,fa->ijab', self.l2, tmp)

//sort!!!  
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


    //tmp   = ndot('me,fiea->mfia', self.x1, self.Hvovv, prefactor=2.0)
    //tmp  -= ndot('me,fiae->mfia', self.x1, self.Hvovv)
    //r_y2 += ndot('mfia,jmbf->ijab', tmp, self.l2)

//sort!!!  
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rpqs, 5, 10, "WAmEf 2(Am,Ef) - (Am,fE) (EA,mf)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 5, 10, 5, 10, 0, "WAmEf 2(Am,Ef) - (Am,fE) (EA,mf)");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z_iajb");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl); 
    global_dpd_->contract244( &X1, &W, &Z, 1, 0, 0, 1, 0); 	 

    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);


//sort!!!
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, prqs, 10, 10, "2 LIjAb - LIjBa (jb,IA)");  //Can we replace?
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "2 LIjAb - LIjBa (jb,IA)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Ztmp");
    global_dpd_->contract444(&Z, &L2, &Z2, 1, 0, 1, 0);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, prqs, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z2);  

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);   

    
    //tmp   = ndot('me,jmna->ejna', self.x1, self.Hooov)
    //r_y2 += ndot('ineb,ejna->ijab', self.l2, tmp)

//  //Here make a test
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "Z_ijkb (tmp)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);    

    global_dpd_->contract244(&X1, &L2, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1); 

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, psrq, 10, 0, "Z(ia,jk)");
    global_dpd_->buf4_close(&Z);

//sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 10, 0, "WMnIe (Me,nI)");  //Can we replace?
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z2(ib,ja)");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "Z(ia,jk)");    
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMnIe (Me,nI)");      
    global_dpd_->contract444(&Z, &W, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, prsq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);


    //tmp   = ndot('me,mjna->ejna', self.x1, self.Hooov)
    //r_y2 += ndot('nieb,ejna->ijab', self.l2, tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "Z_ijkb (tmp)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 10, "Z (ji,kb)");
    global_dpd_->buf4_close(&Z);

//sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, sqrp, 11, 0, "WMnIe (en,IM)");  //Can we replace?
    global_dpd_->buf4_close(&W);

    
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "Z2 (ib,aj)");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "Z (ji,kb)"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 0, 11, 0, 0, "WMnIe (en,IM)"); 
    global_dpd_->contract444(&Z, &W, &Z2, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);   

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, psrq, 0, 5, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z2);
    
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);


    //tmp   = ndot('me,nmba->enba', self.x1, self.l2)
    //r_y2 += ndot('jine,enba->ijab', self.Hooov, tmp)
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "Z (ij,kl)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl); 
    global_dpd_->contract424(&W, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&X1);

    //sort!!! can we reuse?
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
   

    //tmp   = ndot('me,mina->eina', self.x1, self.Hooov, prefactor=2.0)
    //tmp  -= ndot('me,imna->eina', self.x1, self.Hooov)
    //r_y2 -= ndot('eina,njeb->ijab', tmp, self.l2)

//   X1*L2....
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij,ka)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&X1, &L2, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rpqs, 0, 10, "Z (ik,ja)");
    global_dpd_->buf4_close(&Z);

//sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 0, 10, "2WMnIe - WnMIe (MI,nE)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z (jb,ia)");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "Z (ik,ja)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe (MI,nE)");	
    global_dpd_->contract444(&W, &Z, &Z2, 1, 1, 1, 0);   
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, prqs, 0, 5, "Z(ij,ab)");

    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);


    //tmp   = ndot('me,imne->in', self.x1, self.Hooov, prefactor=2.0)
    //tmp  -= ndot('me,mine->in', self.x1, self.Hooov)
    //r_y2 -= ndot('in,jnba->ijab', tmp, self.l2)

//sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 0, 10, "2WMnIe - WnMIe (MI,nE)");    
    global_dpd_->buf4_close(&W);


    global_dpd_->file2_init(&z, PSIF_CC_OEI, irrep, 0, 0, "z_ij");
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe (MI,nE)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->contract422(&W, &X1, &z, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa (jI,bA)"); //I reuse it    
    global_dpd_->contract244(&z, &L2, &Y2new, 1, 0, 0, -1, 1); 
    global_dpd_->file2_close(&z);
    global_dpd_->buf4_close(&L2);  


           Y2_norm = global_dpd_->buf4_dot_self(&Y2new);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\n\tTODAY NORM Y2NEW FINAL!!!!! %20.15f\n", Y2_norm);


    // **** <O|L2(0)|[Hbar(0), X2]|phi^ab_ij> ****

    //Can we reuse this part?
    //global_dpd_->buf4_init(&test, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "test");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, "Z_ijkl");    
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

  
//sort D
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 0, 5, "D (ij|ba)");   
    global_dpd_->buf4_close(&D);  

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, "Z_ijkl");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D (ij|ba)"); //Switch indicies for X?
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&D, &X2, &Z, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&X2);


    //sort!!! can we reuse? Yes! 
    //global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    //global_dpd_->buf4_sort(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa (Ij,bA)");
    //global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa (Ij,bA)");
    global_dpd_->contract444(&Z, &L2, &Y2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);

//-------------------------------------
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &Z, 1, 1, 1, 0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    //sort D
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D (ia|jb)"); 
    global_dpd_->buf4_close(&D);     

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z2 (ib|ja)"); 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ia|jb)");
    global_dpd_->contract444(&Z, &D, &Z2, 0, 0, 1, 0);   
    global_dpd_->buf4_close(&D);

//sort 
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prsq, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

//-------
    //I think we can use the intermediate
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");
    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &Z, 0, 1, 1, 0); 
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

//sort D
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, qrps, 10, 10, "D (ja|ib)"); 
    global_dpd_->buf4_close(&D);     

 
    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z_tmp");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D (ja|ib)"); 
    global_dpd_->contract444(&Z, &D, &Z2, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

//sort 
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prsq, 0, 5, "Z2 (ij,ab)");  
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

//-------------- 
    //I think we can use the intermediate
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_jbia");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)"); //I am reusing 
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    //Lx_ijab
    global_dpd_->contract444(&L2, &X2, &Z, 1, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_ibja_test");  //I am reusing
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LX_jbia");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)"); //I am reusing
    global_dpd_->contract444(&Z, &D, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

   //sort
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, rpsq, 0, 5, "Z2 (ij,ab)"); 
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);
//---------------------

    sprintf(lbl, "G_%s_oo (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&D, &X2, &GMI, 0, 0, 1.0, 0.0); //Can we reuse it??
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "tmp"); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    global_dpd_->contract424(&L2, &GMI, &Z2, 1, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&GMI);

    //sort
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, qpsr, 0, 5, "Z(ij,ab)new");  //Test maybe it is worth it to sort L2?
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z(ij,ab)new");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);


    sprintf(lbl, "G_oo (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&L2, &X2, &GMI, 0, 0, 1.0, 0.0); //Can we reuse it??
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GMI, &Y2new, 1, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);

//
    sprintf(lbl, "G_%s_vv (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&D, &X2, &GAE, 2, 2, -1.0, 0.0); //Can we reuse it??

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract244(&GAE, &L2, &Y2new, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&L2);


//-----------------------
    sprintf(lbl, "G_%s_vv (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");   
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&L2, &X2, &GAE, 2, 2, -1.0, 0.0); //Can we reuse it??
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    //global_dpd_->contract244(&GAE, &D, &Y2new, 1, 3, 1, 1.0, 1.0);
    global_dpd_->contract424(&D, &GAE, &Y2new, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&D);
//-----------------------------------------------------------

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

   //sort
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prsq, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, -1);
    global_dpd_->buf4_close(&Z2);

//---------------------  
    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Z(ib|ja)");
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "Xl (ib|me)"); //I am reusing
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");  
    global_dpd_->contract444(&D, &Z, &Z2, 0, 0, 2, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&D);

    //sort
    global_dpd_->buf4_sort(&Z2, PSIF_CC_LR, prqs, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &Y2new, 1);
    global_dpd_->buf4_close(&Z2);

//    if (params.local)
//        local_filter_T2(&Y2new);
//    else
//        denom2(&Y2new, omega);

    global_dpd_->buf4_close(&Y2new);

//global_dpd_->buf4_close(&Y2);
//global_dpd_->buf4_close(&Y2new);
}

}  // namespace ccresponse
}  // namespace psi
