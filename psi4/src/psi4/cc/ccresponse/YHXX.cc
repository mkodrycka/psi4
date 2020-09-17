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


double Y1HX1X1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, F, z, z1, Z_final;
    dpdbuf4 W, Z, Z2, Y2 ; 
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;
    
    // *** <O|Y1(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep_x, 0, 5, 0, 5, 0, "Z2 (ij|ab)");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");

    global_dpd_->file2_mat_init(&Y1);
    global_dpd_->file2_mat_rd(&Y1);
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z2 is totally symmetric 
        Gij = Gab ^ irrep_x;
        global_dpd_->buf4_mat_irrep_init(&Z2, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z2, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep_x;
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z2.params->rowtot[Gij]; ij++) {
               i = Z2.params->roworb[Gej][ij][0];
               j = Z2.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z2.params->coltot[Gij]; ab++) {
                   a = Z2.params->colorb[Gab][ab][0];
                   b = Z2.params->colorb[Gab][ab][1];
                   Z2.matrix[Gij][ij][ab] -= F.matrix[Gi][i][a] * Y1.matrix[Gj][j][b]; 
                   Z2.matrix[Gij][ij][ab] -= Y1.matrix[Gj][i][a] * F.matrix[Gi][j][b];
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z2, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z2, Gij);
    }

    global_dpd_->file2_mat_close(&Y1);
    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_close(&F);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_scm(&Z, 0);
    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->contract424(&W, &Y1, &Z, 3, 0, 0, -1, 1);

    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe (nM,eI)");   //sort
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (nM,eI)");
    global_dpd_->contract424(&W, &Y1, &Z2, 3, 0, 0, -1, 0);   
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, pqsr, 0, 5, "Z2 (ij,ba)");   //sort
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ba)");
    global_dpd_->buf4_axpy(&Z2, &Z, 1);  
    global_dpd_->buf4_close(&Z2);

//Here Compute out of core!!!!!!!!

/*
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);
*/

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    global_dpd_->contract424(&W, &Y1, &Z2, 1, 1, 1, 1, 0);   

    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");
    global_dpd_->contract244(&Y1, &W, &Z2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&Y1);

    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);


    //-------------------------------------------------------------//

    sprintf(lbl, "Z_%s_Final (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&Z_final, PSIF_CC_OEI, irrep_y, 0, 1, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->dot14(&X1, &Z, &Z_final, 0, 0, 1, 0);

    global_dpd_->buf4_close(&Z);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);

    result = global_dpd_->file2_dot(&Z_final, &X1); 

    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&Z_final);


    return result;
}


double Y2HX1X1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, F, FX, z, z1, Z_final;
    dpdbuf4 W, Z, Z2, Y2 ;
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;

    // *** <O|Y2(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    // Imeplement this part!!!!!!!!!!!


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_scm(&Z, 0);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Y2);  

//outfile->Printf("\n\tResult B1:  %20.15f\n", result);


    return result;
}



double Y2HX2X2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, z_ij, z_ab, F, FX, z, z1; 
    dpdbuf4 X2,Z, D, Y2, Z_final ;
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;

    // *** <O|Y2(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    //tmp = np.einsum("klcd,ijcd->ijkl",self.x2_C,self.y2_A)    
    //tmp = np.einsum("ijkl,ijab->klab",tmp,self.x2_B)
    //self.Bcon1 += 0.5*np.einsum('klab,klab->',tmp,self.Goovv)

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_scm(&Z_final, 0);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z_ijkl");    
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&Y2, &X2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");  
    global_dpd_->contract444(&Z, &D, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);


    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result = 0.5 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);
    

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z (ij,kl)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z, &D, &Z_final, 1, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);


    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);

    result += 0.5 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ijab,ikbd->jkad",self.x2_B,self.y2_A)    
    //tmp = np.einsum("jkad,klcd->jlac",tmp,self.x2_C)
    //self.Bcon1 += np.einsum('jlac,jlac->',tmp,self.Goovv) 

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 1, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij,ab)");   
    global_dpd_->buf4_close(&Z_final); 


    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");    
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    result += global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("klcd,ikdb->licb",self.x2_C,self.y2_A)
    //tmp = np.einsum("licb,ijab->ljca",tmp,self.x2_B)
    //self.Bcon1 += np.einsum('ljca,ljac->',tmp,self.Goovv) 

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z_final (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prsq, 0, 5, "Z_final (ij,ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z_final (ij,ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    result += global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ijab,ijac->bc",self.x2_B,self.Loovv)  
    //tmp = np.einsum("bc,klcd->klbd",tmp,self.x2_C)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ab, 3, 3, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ab, &X2, &Z_final, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ab);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);



    //tmp = np.einsum("ijab,ikab->jk",self.x2_B,self.Loovv)  
    //tmp = np.einsum("jk,klcd->jlcd",tmp,self.x2_C)
    //self.Bcon1 -= np.einsum("jlcd,jlcd->",tmp,self.y2_A)

    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ij, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ij, &X2, &Z_final, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ij);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ikbc,klcd->ilbd",self.Loovv,self.x2_C)
    //tmp = np.einsum("ilbd,ijab->jlad",tmp,self.x2_B)
    //self.Bcon1 -= np.einsum("jlad,jlad->",tmp,self.y2_A)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&D, &X2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, rpsq, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&Y2);



    //tmp = np.einsum("ijab,jlbc->ilac",self.x2_B,self.y2_A)
    //tmp = np.einsum("ilac,klcd->ikad",tmp,self.x2_C)
    //self.Bcon1 -= np.einsum("ikad,ikad->",tmp,self.Loovv)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);


    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");

    result -= global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&D);


    //tmp = np.einsum("klca,klcd->ad",self.Loovv,self.x2_C)
    //tmp = np.einsum("ad,ijdb->ijab",tmp,self.y2_A)
    //self.Bcon1 -= np.einsum("ijab,ijab->",tmp,self.x2_B)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&D, &X2, &z_ab, 3, 3, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ab, &Y2, &Z_final, 1, 2, 1, 1, 0);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->file2_close(&z_ab);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("kicd,klcd->il",self.Loovv,self.x2_C)
    //tmp = np.einsum("ijab,il->ljab",self.x2_B,tmp)
    //self.Bcon1 -= np.einsum("ljab,ljab->",tmp,self.y2_A)

    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ij, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ij, &X2, &Z_final, 1, 0, 0, 1, 0);

    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ij);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    
    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("klcd,ikac->lida",self.x2_C,self.y2_A)
    //tmp = np.einsum("lida,jlbd->ijab",tmp,self.Loovv)
    //self.Bcon1 += 2.*np.einsum("ijab,ijab->",tmp,self.x2_B)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&Z, &D, &Z_final, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result += 2 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&X2);


    outfile->Printf("\n\tResult B1:  %20.15f\n", result);


    return result;
}


double YHXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		     const char *pert_z, int irrep_z, double omega_z) {

    double hyper = 0.0;

 
    // *** <O|Y1(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    hyper += Y1HX1X1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 

   // ***  <O|Y2(A)|[[Hbar(0),X1(B)],X1(C)]|0> ***

    hyper += Y2HX1X1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

   // ***  <O|Y2(A)|[[Hbar(0),X2(B)],X2(C)]|0> ***

    hyper += Y2HX2X2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);    

       

    //outfile->Printf("\n\tHYPER B1:  %20.15f\n", hyper);

    return hyper;
}

}  // namespace ccresponse
}  // namespace psi
