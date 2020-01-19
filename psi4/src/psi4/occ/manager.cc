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

#include "occwave.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

#include <cmath>

using namespace psi;

namespace psi {
namespace occwave {

//======================================================================
//             OMP2 Manager
//======================================================================
void OCCWave::omp2_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("MP2 Energy");
    omp2_mp2_energy();
    timer_off("MP2 Energy");
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;
    if (ip_poles == "TRUE") omp2_ip_poles();
    if (ep_ip_poles == "TRUE") ep2_ip();

    mp2_postprocessing();

    omp2_response_pdms();
    gfock();
    idp();
    mograd();
    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOEscsmp2 - EscfoT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard MP2 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        set_t2_amplitudes_mp2();
        conver = 1;
        if (dertype == "FIRST") {
            omp2_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        omp2_mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;

        // S2
        // if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();

        // Green's function
        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
            }
            omp2_ip_poles();
        }

        if (ep_ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
            }
            ep2_ip();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp2_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP2 FINAL RESULTS ========================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tSCS-OMP2 Total Energy (a.u.)       : %20.14f\n", Escsmp2);
        outfile->Printf("\tSOS-OMP2 Total Energy (a.u.)       : %20.14f\n", Esosmp2);
        outfile->Printf("\tSCSN-OMP2 Total Energy (a.u.)      : %20.14f\n", Escsnmp2);
        outfile->Printf("\tSCS-OMP2-VDW Total Energy (a.u.)   : %20.14f\n", Escsmp2vdw);
        outfile->Printf("\tSOS-PI-OMP2 Total Energy (a.u.)    : %20.14f\n", Esospimp2);
        outfile->Printf("\tOMP2 Correlation Energy (a.u.)     : %20.14f\n", Emp2L - Escf);
        outfile->Printf("\tEomp2 - Eref (a.u.)                : %20.14f\n", Emp2L - Eref);
        outfile->Printf("\tOMP2 Total Energy (a.u.)           : %20.14f\n", Emp2L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP2 TOTAL ENERGY"] = Emp2L;
        Process::environment.globals["OMP2 TOTAL ENERGY"] = Emp2L;
        Process::environment.globals["SCS-OMP2 TOTAL ENERGY"] = Escsmp2;
        Process::environment.globals["SOS-OMP2 TOTAL ENERGY"] = Esosmp2;
        Process::environment.globals["SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
        Process::environment.globals["SCS-OMP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
        Process::environment.globals["SOS-PI-OMP2 TOTAL ENERGY"] = Esospimp2;

        // LAB: variables_ and energy_ here are what I vouch for and test.
        //      The P::e.globals will diminish and go into the West and be
        //      replaced by qcvar formulas computed py-side from wfn vars.
        variables_["CURRENT REFERENCE ENERGY"] = Escf;

        variables_["OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["SCS-OMP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-OMP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-OMP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-OMP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        // if scs on
        if (do_scs == "TRUE") {
            if (scs_type_ == "SCS") {
                variables_["CURRENT ENERGY"] = Escsmp2;
            }

            else if (scs_type_ == "SCSN") {
                variables_["CURRENT ENERGY"] = Escsnmp2;
            }

            else if (scs_type_ == "SCSVDW") {
                variables_["CURRENT ENERGY"] = Escsmp2vdw;
            }
        }

        // else if sos on
        else if (do_sos == "TRUE") {
            if (sos_type_ == "SOS") {
                variables_["CURRENT ENERGY"] = Esosmp2;
            }

            else if (sos_type_ == "SOSPI") {
                variables_["CURRENT ENERGY"] = Esospimp2;
            }
        } else {
            variables_["CURRENT ENERGY"] = Emp2L;
        }
        energy_ = variables_["CURRENT ENERGY"];
        variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp2_manager

//======================================================================
//             MP2 Manager
//======================================================================
void OCCWave::mp2_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
    } else {
        if (reference_ == "RESTRICTED")
            trans_ints_rmp2();
        else if (reference_ == "UNRESTRICTED" && reference == "ROHF")
            trans_ints_uhf();
        else if (reference_ == "UNRESTRICTED" && reference != "ROHF")
            trans_ints_ump2();
    }
    timer_off("trans_ints");
    // ROHF REF
    if (reference == "ROHF") {
        timer_on("T1(1)");
        t1_1st_sc();
        timer_off("T1(1)");
    }  // end if (reference == "ROHF")
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    Eref = Escf;
    timer_on("MP2 Energy");
    omp2_mp2_energy();
    timer_off("MP2 Energy");
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;
    if (ip_poles == "TRUE") omp2_ip_poles();
    if (ep_ip_poles == "TRUE") ep2_ip();

    mp2_postprocessing(reference == "ROHF");

    // Why is the below line commented?
    // if (reference == "ROHF") Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    if (do_scs == "TRUE") {
        if (scs_type_ == "SCS") {
            variables_["CURRENT ENERGY"] = Escsmp2;
        }
        else if (scs_type_ == "SCSN") {
            variables_["CURRENT ENERGY"] = Escsnmp2;
        }
        else if (scs_type_ == "SCSVDW") {
            variables_["CURRENT ENERGY"] = Escsmp2vdw;
        }
    }
    else if (do_sos == "TRUE") {
        if (sos_type_ == "SOS") {
            variables_["CURRENT ENERGY"] = Esosmp2;
        }
        else if (sos_type_ == "SOSPI") {
            variables_["CURRENT ENERGY"] = Esospimp2;
        }
    } else {
        variables_["CURRENT ENERGY"] = Emp2;
    }
    energy_ = variables_["CURRENT ENERGY"];
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];


    /* updates the wavefunction for checkpointing */
    name_ = "MP2";

    // S2
    // if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp2_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
            // outfile->Printf("\tAn EKT computation for a non-OO method requested. Analytic gradients will not be
            // computed! \n");
            // tstop();
            // exit(EXIT_SUCCESS);
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end mp2_manager

//======================================================================
//             OMP3 Manager
//======================================================================
void OCCWave::omp3_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    omp3_mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp3_postprocessing();

    omp3_response_pdms();
    gfock();
    idp();
    mograd();
    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard MP3 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        set_t2_amplitudes_mp2();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            omp3_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        omp3_mp2_energy();
        mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
                t2_2nd_sc();
            }
            omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp3_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();

        mp3_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP3 FINAL RESULTS ========================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tSCS-OMP3 Total Energy (a.u.)       : %20.14f\n", Escsmp3);
        outfile->Printf("\tSOS-OMP3 Total Energy (a.u.)       : %20.14f\n", Esosmp3);
        outfile->Printf("\tSCSN-OMP3 Total Energy (a.u.)      : %20.14f\n", Escsnmp3);
        outfile->Printf("\tSCS-OMP3-VDW Total Energy (a.u.    : %20.14f\n", Escsmp3vdw);
        outfile->Printf("\tSOS-PI-OMP3 Total Energy (a.u.)    : %20.14f\n", Esospimp3);
        outfile->Printf("\tOMP3 Correlation Energy (a.u.)     : %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEomp3 - Eref (a.u.)                : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tOMP3 Total Energy (a.u.)           : %20.14f\n", Emp3L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP3 TOTAL ENERGY"] = Emp3L;
        Process::environment.globals["SCS-OMP3 TOTAL ENERGY"] = Escsmp3;
        Process::environment.globals["SOS-OMP3 TOTAL ENERGY"] = Esosmp3;
        Process::environment.globals["SCSN-OMP3 TOTAL ENERGY"] = Escsnmp3;
        Process::environment.globals["SCS-OMP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
        Process::environment.globals["SOS-PI-OMP3 TOTAL ENERGY"] = Esospimp3;

        variables_["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;
        Process::environment.globals["SCS-OMP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
        Process::environment.globals["SOS-OMP3 CORRELATION ENERGY"] = Esosmp3 - Escf;
        Process::environment.globals["SCSN-OMP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
        Process::environment.globals["SCS-OMP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
        Process::environment.globals["SOS-PI-OMP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        if (do_scs == "TRUE") {
            if (scs_type_ == "SCS") {
                variables_["CURRENT ENERGY"] = Escsmp3;
            }

            else if (scs_type_ == "SCSN") {
                variables_["CURRENT ENERGY"] = Escsnmp3;
            }

            else if (scs_type_ == "SCSVDW") {
                variables_["CURRENT ENERGY"] = Escsmp3vdw;
            }
        }
        else if (do_sos == "TRUE") {
            if (sos_type_ == "SOS") {
                variables_["CURRENT ENERGY"] = Esosmp3;
            }

            else if (sos_type_ == "SOSPI") {
                variables_["CURRENT ENERGY"] = Esospimp3;
            }
        } else {
            variables_["CURRENT ENERGY"] = Emp3;
        }
        variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp3_manager

//======================================================================
//             MP3 Manager
//======================================================================
void OCCWave::mp3_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    omp3_mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp3_postprocessing();

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    if (do_scs == "TRUE") {
        if (scs_type_ == "SCS") {
            variables_["CURRENT ENERGY"] = Escsmp3;
        }

        else if (scs_type_ == "SCSN") {
            variables_["CURRENT ENERGY"] = Escsnmp3;
        }

        else if (scs_type_ == "SCSVDW") {
            variables_["CURRENT ENERGY"] = Escsmp3vdw;
        }
    }
    else if (do_sos == "TRUE") {
        if (sos_type_ == "SOS") {
            variables_["CURRENT ENERGY"] = Esosmp3;
        }

        else if (sos_type_ == "SOSPI") {
            variables_["CURRENT ENERGY"] = Esospimp3;
        }
    } else {
        variables_["CURRENT ENERGY"] = Emp3;
    }
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp3_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end mp3_manager

//======================================================================
//             OCEPA Manager
//======================================================================
void OCCWave::ocepa_manager() {
    mo_optimized = 0;  // means MOs are not optimized yet.
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    ocepa_mp2_energy();
    timer_off("MP2 Energy");
    Ecepa = Emp2;
    EcepaL = Ecepa;
    EcorrL = Ecorr;
    EcepaL_old = Ecepa;

    mp2_postprocessing();

    ocepa_response_pdms();
    gfock();
    idp();
    mograd();
    if (rms_wog > tol_grad)
        occ_iterations();
    else {
        orbs_already_opt = 1;
        outfile->Printf("\n\tOrbitals are already optimized, switching to the canonical CEPA computation... \n");

        cepa_iterations();
    }

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard CEPA computation... \n");

        ref_energy();
        cepa_energy();
        Ecepa_old = EcepaL;
        cepa_iterations();
        if (dertype == "FIRST") {
            ocepa_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        cepa_energy();
        if (orbs_already_opt == 1) EcepaL = Ecepa;

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_opt == 1) {
                ocepa_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OCEPA FINAL RESULTS ========================================= \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tSCS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Escscepa);
        outfile->Printf("\tSOS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Esoscepa);
        outfile->Printf("\tOCEPA(0) Correlation Energy (a.u.) : %20.14f\n", EcepaL - Escf);
        outfile->Printf("\tEocepa - Eref (a.u.)               : %20.14f\n", EcepaL - Eref);
        outfile->Printf("\tOCEPA(0) Total Energy (a.u.)       : %20.14f\n", EcepaL);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OLCCD TOTAL ENERGY"] = EcepaL;
        Process::environment.globals["SCS-OLCCD TOTAL ENERGY"] = Escscepa;
        Process::environment.globals["SOS-OLCCD TOTAL ENERGY"] = Esoscepa;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;

        variables_["OLCCD CORRELATION ENERGY"] = EcepaL - Escf;
        Process::environment.globals["SCS-OLCCD CORRELATION ENERGY"] = Escscepa - Escf;
        Process::environment.globals["SOS-OLCCD CORRELATION ENERGY"] = Esoscepa - Escf;

        if (do_scs == "TRUE") {
            variables_["CURRENT ENERGY"] = Escscepa;
        }
        else if (do_sos == "TRUE") {
            variables_["CURRENT ENERGY"] = Esoscepa;
        } else {
            variables_["CURRENT ENERGY"] = EcepaL;
        }
        variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end ocepa_manager

//======================================================================
//             CEPA Manager
//======================================================================
void OCCWave::cepa_manager() {
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    ocepa_mp2_energy();
    timer_off("MP2 Energy");
    Ecepa = Emp2;
    Ecepa_old = Emp2;

    mp2_postprocessing();

    // Perform CEPA iterations
    cepa_iterations();

    outfile->Printf("\n");
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\t================ CEPA FINAL RESULTS ========================================== \n");
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tCEPA(0) Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCEPA(0) Total Energy (a.u.)        : %20.14f\n", Ecepa);
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");

    // Set the global variables with the energies
    variables_["LCCD TOTAL ENERGY"] = Ecepa;
    variables_["LCCD CORRELATION ENERGY"] = Ecorr;
    variables_["CURRENT ENERGY"] = Ecepa;
    variables_["CURRENT REFERENCE ENERGY"] = Eref;
    variables_["CURRENT CORRELATION ENERGY"] = Ecorr;
    // EcepaL = Ecepa;

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        ocepa_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }
}  // end cepa_manager

//======================================================================
//             OMP2.5 Manager
//======================================================================
void OCCWave::omp2_5_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    omp3_mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp2p5_postprocessing();

    omp3_response_pdms();
    gfock();
    // ccl_energy();
    idp();
    mograd();
    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard MP2.5 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        set_t2_amplitudes_mp2();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            omp3_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        omp3_mp2_energy();
        mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
                t2_2nd_sc();
            }
            omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp3_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();
        mp2p5_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP2.5 FINAL RESULTS ======================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tOMP2.5 Correlation Energy (a.u.)   : %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEomp2.5 - Eref (a.u.)              : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tOMP2.5 Total Energy (a.u.)         : %20.14f\n", Emp3L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP2.5 TOTAL ENERGY"] = Emp3L;
        variables_["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["CURRENT ENERGY"] = Emp3L;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp2.5_manager

//======================================================================
//             MP2.5 Manager
//======================================================================
void OCCWave::mp2_5_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    omp3_mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp2p5_postprocessing();
    variables_["CURRENT ENERGY"] = Emp3L;
    variables_["CURRENT REFERENCE ENERGY"] = Eref;
    variables_["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp3_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end omp2.5_manager
}
}  // End Namespaces
