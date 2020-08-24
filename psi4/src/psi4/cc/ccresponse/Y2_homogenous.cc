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

void Y2_homogenous_build(const char *pert, int irrep, double omega) {
    dpdfile2 Y1, z, F, GAE, GMI, t1;
    dpdbuf4 Y2, Y2new, Y2inhomo, W, D,  Z, Z1, Z2, B;
    char lbl[32];
    int Gej, Gab, Gij, Ge, Gj, Gi, nrows, length, E, e, II;
    int Gbm, Gfe, bm, b, m, Gb, Gm, Gf, M, fe, f, ef, ncols;
    double *Y;
    dpdbuf4 S, A, B_s;
    int ij, Gc, C, c, cc;
    int rows_per_bucket, nbuckets, row_start, rows_left, nlinks;
    psio_address next;
    double **Y_diag, **B_diag;
    double Y2_norm;

    //Set of homogenous terms 

    //a factor of 0.5 because teh amplitudes are symmetric
   
    sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega); 
    global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 

    global_dpd_->buf4_scm(&Y2new, 0); //Do I need this??  

    //Add Inhomogenous terms
    sprintf(lbl, "Inhomo Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2inhomo, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 

    global_dpd_->buf4_axpy(&Y2inhomo, &Y2new, 1);
    global_dpd_->buf4_close(&Y2inhomo);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_axpy(&Y2, &Y2new, 0.5*omega);    //Make sure about 0.5

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract244(&F, &Y2, &Y2new, 0, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F); 


    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract244(&F, &Y2, &Y2new, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&F);


    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&W, &Y2, &Y2new, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&W);     
	

    //Y2*Hvvvv COmpute out of core!!! It is not working!!!!!
    /* RHS += Wefab*Yijef  */
    //This part of the code must be replaced! See WefabL2.cc  
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "ZAbIj");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract444(&B, &Y2, &Z, 0, 0, 0.5, 0);
    global_dpd_->buf4_close(&B);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rspq, 0, 5, "ZIjAb_W");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "ZIjAb_W"); 
    //global_dpd_->buf4_axpy(&Z, &Y2new, 1);   
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); //Make it out of core
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl); 
    global_dpd_->contract244(&Y1, &W, &Y2new, 1, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (nM,eI)"); 
    global_dpd_->contract424(&W, &Y1, &Y2new, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Y1);
    global_dpd_->buf4_close(&W);


    //It is NOT WORKING!!!!!!!!!!!!
    //Build 2WMbEj - WMbjE
/*
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WMbEj 2 W(Mb,Ej) - W(Mb,jE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 10, 11, "WMbEj 2 W(Mb,Ej) - W(Mb,jE)", -1);
    global_dpd_->buf4_close(&W);

outfile->Printf("\tNorm of Here1\n");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj 2 W(Mb,Ej) - W(Mb,jE)");
outfile->Printf("\tNorm of Here2\n");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 10, 10, "WMbEj 2 W(Mb,jE) - W(Mb,jE)");
outfile->Printf("\tNorm of Here3\n");    

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj 2 W(Mb,jE) - W(Mb,jE)");
outfile->Printf("\tNorm of Here4\n");
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Ziajb"); 
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 10, 10, 10, 10, lbl);     
outfile->Printf("\tNorm of Here5\n");
    global_dpd_->contract444(&W, &Y2, &Z, 0, 1.0, 1.0, 0);
outfile->Printf("\tNorm of Here6\n");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Y2);
outfile->Printf("\tNorm of Here7\n");   
 
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Ziajb"); 
outfile->Printf("\tNorm of Here8\n");
    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 5, "Zijab_sorted");
    global_dpd_->buf4_close(&Z);
outfile->Printf("\tNorm of Here9\n");  

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Zijab_sorted");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1.0);  
outfile->Printf("\tNorm of Here10\n");
    global_dpd_->buf4_close(&Z);
*/

    //(2*Hovov-Hovvo)*Y2 Thi is not working!!!!!!
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    //global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 0, 5, "2 W(Mj,Eb) + W(MJ,Eb)");   
    //global_dpd_->buf4_close(&W);

/*
//It is not working!!!!!!!!!!!!!!
   // r_y2 -= ndot('mieb,jeam->ijab', self.y2, self.Hovvo)
   // sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 10, 11, "WMEbj");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Ziajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMEbj");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Y2, &W, &Z, 0, 1, 1.0, 0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 5, "Zijab_sorted");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Zijab_sorted");

           Y2_norm = global_dpd_->buf4_dot_self(&Z);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\tNorm of Zijab_sorted start %20.15f\n", Y2_norm);


    global_dpd_->buf4_axpy(&Z, &Y2new, -1.0);
    global_dpd_->buf4_close(&Z);
*/
//------------------------------------------------------------------
    sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl); 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract244(&GAE, &D, &Y2new, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);


    sprintf(lbl, "G_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract244(&GMI, &D, &Y2new, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);    

    if (params.local)
        local_filter_T2(&Y2new);
    else
        denom2(&Y2new, omega);
    global_dpd_->buf4_close(&Y2new);

//global_dpd_->buf4_close(&Y2);
//global_dpd_->buf4_close(&Y2new);
}

}  // namespace ccresponse
}  // namespace psi
