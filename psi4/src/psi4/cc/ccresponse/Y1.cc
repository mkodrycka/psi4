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

   Computes inhomogenous terms appearing in Y1 equations.
   
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

void local_filter_T1(dpdfile2 *);
void lambda_residuals();
void L2HX2(const char *pert, int irrep, double omega);

void Y1_inhomogenous_build(const char *pert, int irrep, double omega) {
    dpdfile2 ME, newYIA, newYia, YIA, Yia;
    dpdfile2 GAE, GMI;
    int Gim, Gi, Gm, Ga, Gam, nrows, ncols, A, a, am;
    int Gei, ei, e, i, Gef, Ge, Gf, E, I, af, fa, f;
    int GW, GX1, GZ, Gej, Gab, Gij, Gj;
    int num_j, num_i, num_e, nlinks;
    double *X;
    dpdfile2 F, z1, z2;
    dpdbuf4 W, WL, D, X2, Z2, Z3, lx_iajb, X2test, L2test, LIjAb; 
    dpdfile2 Y1, Y1new, mu1, L1, lt, lx, lx_AB, X1;
    dpdbuf4 L2, Z, mu2, Hx_ijab, lx_ijab;
    char lbl[32];
    double Y1_norm;
    double *Y;
    dpdfile2 test;

    sprintf(lbl, "Inhomo Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

    /*** Mu * L1 ***/ 
    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_axpy(&mu1, &Y1new, 2, 0);  
    global_dpd_->file2_close(&mu1);

    /*** L1 * MuBAR + L2 * MuBAR ***/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    sprintf(lbl, "%sBAR_MI", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract222(&mu1, &L1, &Y1new, 0, 1, -1, 1.0);
    global_dpd_->file2_close(&mu1);

    sprintf(lbl, "%sBAR_AE", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->contract222(&L1, &mu1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&mu1);
    global_dpd_->file2_close(&L1);

    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&lt, PSIF_CC_OEI, 0, 0, 0, "Lt_IJ");
    global_dpd_->contract222(&lt, &mu1, &Y1new, 0, 1, 2, 1.0);  

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    sprintf(lbl, "%sBAR_MbIj", pert, omega);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "%sBAR_MbIj", pert);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);

    // <O|[Hbar(0), X1]|0>
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &Y1new, 0, 0, 2, 1); 
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&X1);   

    // <O|L1(0)|[Hbar(0), X1]|0>
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->contract424(&W, &L1, &Z, 3, 0, 0, -1, 0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 0, 5, "Z (ji,ba)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij,ab)");   
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ji,ba)");  

    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->buf4_init(&WL, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WL_p1 (ij,ab)");
    global_dpd_->buf4_axpy(&WL, &Z, 1);
    global_dpd_->buf4_close(&WL);
  
    global_dpd_->buf4_init(&WL, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WL_p2 (ij,ab)");
    global_dpd_->buf4_axpy(&WL, &Z, 1);
    global_dpd_->buf4_close(&WL);    
 
    global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&L1); 
    global_dpd_->buf4_close(&Z);

    // <O|L2(0)|[Hbar(0), X1]|0>
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract444(&W, &L2, &Z, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ji,ab) - (ji,ba)"); 
    global_dpd_->buf4_axpy(&WL, &Z, 0.5);
    global_dpd_->buf4_close(&WL); 

    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ba) - (ij,ab)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_axpy(&WL, &Z, 0.5); 
    global_dpd_->dot23(&X1, &Z, &Y1new, 0, 0, 1, 1);	
    global_dpd_->buf4_close(&WL);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "GAED (ij,ab)");
    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, 0, 1, 1, "GAE");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GAE, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&D);    
      
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->dot24(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->dot13(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "GMID (ij,ab)");
    
    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, 0, 0, 0, "GMI");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GMI, &Z, 1, 0, 1, 1, 0);
    global_dpd_->dot13(&X1, &Z, &Y1new, 0, 0, -1, 1);
    global_dpd_->dot24(&X1, &Z, &Y1new, 0, 0, -1, 1);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);    

    // Type-II L2 residual 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual II");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&L2, &X1, &Y1new, 0, 0, -1, 1);  
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual IItest");
    global_dpd_->contract422(&L2, &X1, &Y1new, 0, 0, -1, 1); 
    global_dpd_->file2_close(&X1);

    //# <O|L1(0)|[Hbar(0), X2]|phi^a_i>
    sprintf(lbl, "Z_%s_ME", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->dot24(&L1, &X2, &z1, 0, 0, 2, 0); 
    
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);  
    sprintf(lbl, "Z2_%s_ME", pert);
    global_dpd_->file2_init(&z2, PSIF_CC_TMP0, irrep, 0, 1, lbl);
    global_dpd_->dot23(&L1, &X2, &z2, 0, 0, -1, 0);
    global_dpd_->file2_axpy(&z2, &z1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z2);
 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot24(&z1, &D, &Y1new, 0, 0, 1, 1); 
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_MN", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 0, 0, 1.0, 0.0);
    global_dpd_->contract222(&z1, &L1, &Y1new, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_AE", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 1, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 2, 2, -1.0, 0.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);
    global_dpd_->contract222(&L1, &z1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&L1);

    // <O|L2(0)|[Hbar(0), X2]|0> 
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 0, 0, "Lx_IJ");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&L2, &X2, &lx, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 0, 0, "Lx_IJ");    
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME"); 
    global_dpd_->contract222(&lx, &F, &Y1new, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);  

    // Lijab * Xijab -> Lx_AB 
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 1, 1, "Lx_AB");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&X2, &L2, &lx, 2, 2, 1.0, 0.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 1, 1, "Lx_AB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");	
    global_dpd_->contract222(&F, &lx, &Y1new, 0, 1, -1.0, 1.0);   
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)"); 
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);
 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (am,fe)"); //Compute this part out of core
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);  
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);    
    global_dpd_->buf4_close(&W);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&Hx_ijab, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "Hx_ijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
   
    global_dpd_->contract444(&X2, &W, &Hx_ijab, 0, 0, 1, 0);
    global_dpd_->contract442(&Hx_ijab, &L2, &Y1new, 3, 3, -1, 1);

    global_dpd_->buf4_close(&Hx_ijab);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl); 
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, "LX (ia,jb)");

    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");	
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 1, 1, "Lx_AB");
    global_dpd_->dot13(&lx, &W, &Y1new, 1, 0, 1.0, 1.0); 
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&lx);    

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, irrep, 0, 0, 0, 0, 0, "Lx_ijkl");
    
    global_dpd_->contract444(&L2, &X2, &lx_ijab, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract442(&lx_ijab, &W, &Y1new, 1, 3, 1, 1);
    global_dpd_->buf4_close(&lx_ijab);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");

    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);

    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_3");  
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");  

    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (nM,eI)");
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 0, 0, "Lx_IJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->dot14(&lx, &W, &Y1new, 1, 0, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&W);
   
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 1, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_iajb);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)"); 
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, -1, 1);   
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&Y1new);

    return;
}

}  // namespace ccresponse
}  // namespace psi
