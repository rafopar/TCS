/*
 * File:   AnaPiPiProt.cc
 * Author: rafopar
 *
 * Analyzes a (skimmed) hipo file, e.g. OutSkims/005163/PionPairSkim_5163_All.hipo,
 * and selects events with exactly one pi-, one pi+ and one proton.
 *
 * Particle selection, RecParticle/RecEvent usage and the vt computation follow
 * ~/work/git/DDVCS/elDDVCS/AnaElDDVCS.cc (proton logic at lines 135-163):
 *   - proton: chi2pid in (protChi2PIDMin, protChi2PIDMax); vt uses the highest
 *     priority valid scintillator response (FTOF1B > FTOF1A > FTOF2 > CTOF) and
 *     beta computed from the proton mass.
 *   - pi-/pi+: forward-detector status cut (abs(status) in [2000,4000)); vt is
 *     computed the same way as for the proton, but with beta computed from the
 *     charged pion mass.
 *
 * Usage: AnaPiPiProt.exe <input.hipo>
 *
 * Output histograms are written to Hists/Hists_<input basename>.root, unless the
 * input basename matches the MC naming convention MC_Rho_pipi_<jobID>_<index>.hipo
 * (e.g. MC_Rho_pipi_11833_11.hipo), in which case output goes to
 * Hists/MC/Job_<jobID>/Hists_PionPair_<index>.root instead.
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <regex>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Math/Vector4D.h>
#include <TMath.h>
#include <TSystem.h>

// ===== Hipo headers =====
#include <reader.h>
#include <dictionary.h>
#include <RecParticle.h>
#include <RecEvent.h>
#include <MCEvent.h>

using namespace std;

namespace {

    const double vlight = 29.979; // cm/ns
    const double Mpi = 0.13957;
    const double Mprot = 0.9383;

    const double Eb = 10.6; // GeV, beam energy

    const double Mx2_Max = 0.5;
    const double Mx2_Min = -0.5;

    const double Q2_Max = 0.15;

    const double emFT_vt_max = 2.;
    const double emFT_vt_min = -2.;

    const double pim_vt_max = 0.8;
    const double pim_vt_min = -0.8;

    const double pip_vt_max = 0.8;
    const double pip_vt_min = -0.8;

    const double pim_vz_max = 3.;
    const double pim_vz_min = -8.;

    const double pip_vz_max = 2.;
    const double pip_vz_min = -10.;

    const double prot_vz_Max = 1.;
    const double prot_vz_Min = -7;

    const int STATUS_MIN = 2000;
    const int STATUS_MAX = 4000;

    const double protChi2PIDMin = -10.;
    const double protChi2PIDMax = 10.;

    // Highest-priority valid scintillator response for a particle, as in
    // AnaElDDVCS.cc's proton vt logic (FTOF1B > FTOF1A > FTOF2 > CTOF).
    const ScintillatorResponse* pickScintillator(const RecPart& p) {
        if (p.FTOF1B().valid) return &p.FTOF1B();
        if (p.FTOF1A().valid) return &p.FTOF1A();
        if (p.FTOF2().valid) return &p.FTOF2();
        if (p.CTOF().valid) return &p.CTOF();
        return nullptr;
    }

    // vt assuming particle mass `mass`, as in AnaElDDVCS.cc lines 155-157.
    bool computeVt(const RecPart& p, double mass, double& vt) {
        const ScintillatorResponse* sc = pickScintillator(p);
        if (sc == nullptr) {
            return false;
        }

        double beta = p.p() / sqrt(p.p() * p.p() + mass * mass);
        vt = sc->time - sc->path / (beta * vlight) - p.vt();
        return true;
    }

    // User provides vz, so that the path will be calculated as the distance between the cluster (x,y,z) and (0., 0., vz)
    bool computeVtFT(const RecPart &p, double vz, double& vt) {

        if ( !p.FTCal().valid ) {
            return false;
        }

        double path = sqrt(p.FTCal().x*p.FTCal().x + p.FTCal().y*p.FTCal().y + (p.FTCal().z - vz)*(p.FTCal().z - vz) );

        vt  = p.FTCal().time - path/vlight;
        return true;
    }

    // Invariant mass of the generated pi+/pi- pair (MC truth), from the first
    // pid==211 and pid==-211 entries in MC::Particle. Returns false (M_pipi
    // left untouched) if either isn't found -- e.g. real data, which has no
    // MC::Particle rows.
    bool computeMC_Minv_pippim(const MCEvent& mcEvent, double& M_pipi) {
        const MCPart* mcPip = nullptr;
        const MCPart* mcPim = nullptr;
        for (const auto& p : mcEvent.Particles()) {
            if (p.pid() == 211 && mcPip == nullptr) {
                mcPip = &p;
            } else if (p.pid() == -211 && mcPim == nullptr) {
                mcPim = &p;
            }
        }
        if (mcPip == nullptr || mcPim == nullptr) {
            return false;
        }

        ROOT::Math::PxPyPzEVector L_mcPip, L_mcPim;
        L_mcPip.SetPxPyPzE(mcPip->px(), mcPip->py(), mcPip->pz(),
                sqrt(mcPip->p() * mcPip->p() + mcPip->mass() * mcPip->mass()));
        L_mcPim.SetPxPyPzE(mcPim->px(), mcPim->py(), mcPim->pz(),
                sqrt(mcPim->p() * mcPim->p() + mcPim->mass() * mcPim->mass()));

        M_pipi = (L_mcPip + L_mcPim).M();
        return true;
    }
}

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input.hipo>" << endl;
        cout << "Selects events with exactly one pi-, one pi+ and one proton." << endl;
        exit(1);
    }

    const std::string inputFile = argv[1];

    std::string base = inputFile;
    size_t slashPos = base.find_last_of('/');
    if (slashPos != std::string::npos) {
        base = base.substr(slashPos + 1);
    }
    size_t dotPos = base.find_last_of('.');
    std::string stem = (dotPos != std::string::npos) ? base.substr(0, dotPos) : base;

    // MC files are named MC_Rho_pipi_<jobID>_<index>.hipo, e.g.
    // MC_Rho_pipi_11833_11.hipo -> jobID=11833, index=11. These get their own
    // output layout, separate from the flat Hists/Hists_<stem>.root used for data.
    static const std::regex mcNameRe("^MC_Rho_pipi_(\\d+)_(\\d+)\\.hipo$");
    std::smatch mcMatch;

    std::string outputFile;
    if (std::regex_match(base, mcMatch, mcNameRe)) {
        const std::string& jobID = mcMatch[1];
        const std::string& index = mcMatch[2];
        const std::string outDir = "Hists/MC/Job_" + jobID;
        gSystem->Exec(("mkdir -p " + outDir).c_str());
        outputFile = outDir + "/Hists_PionPair_" + index + ".root";
    } else {
        gSystem->Exec("mkdir -p Hists");
        outputFile = "Hists/Hists_" + stem + ".root";
    }

    TFile file_out(outputFile.c_str(), "Recreate");

    TH1D h_vz_pip("h_vz_pip", "", 200, -25., 25.);
    TH1D h_vz_pim("h_vz_pim", "", 200, -25., 25.);
    TH1D h_vz_prot("h_vz_prot", "", 200, -25., 25.);

    TH1D h_chi2PID_pip("h_chi2PID_pip", "", 400, -20., 20.);
    TH1D h_chi2PID_pim("h_chi2PID_pim", "", 400, -20., 20.);
    TH1D h_chi2PID_prot("h_chi2PID_prot", "", 400, -20., 20.);

    TH1D h_chi2PID_pip_FS("h_chi2PID_pip_FS", "", 400, -20., 20.); // When the Final State (FS) is selected
    TH1D h_chi2PID_pim_FS("h_chi2PID_pim_FS", "", 400, -20., 20.); // When the Final State (FS) is selected
    TH1D h_chi2PID_prot_FS("h_chi2PID_prot_FS", "", 400, -20., 20.); // When the Final State (FS) is selected

    TH1D h_vt_pip("h_vt_pip", "", 400, -10., 10.);
    TH1D h_vt_pim("h_vt_pim", "", 400, -10., 10.);
    TH1D h_vt_prot("h_vt_prot", "", 400, -10., 10.);
    TH1D h_vtem_FT("h_vtem_FT", "", 400, -10., 10.);

    TH1D h_Minv_pippim("h_Minv_pippim", "", 200, 0., 3.);
    TH1D h_Minv_pippim_Mx2Cut("h_Minv_pippim_Mx2Cut", "", 200, 0., 3.);

    // MC truth (generated pi+/pi- pair) invariant mass, at three points of the
    // reconstructed-event selection: every event regardless of whether any
    // RecPart was found, once exactly one reconstructed pi-/pi+/proton is
    // selected, and once that selection also passes the Mx2 cut. Filled from
    // the same MC truth value each time -- see computeMC_Minv_pippim -- so
    // these histograms show reconstruction/selection efficiency vs. M(pi+pi-),
    // not three different mass calculations. Stay empty (never filled) for
    // real data, which has no MC::Particle rows.
    TH1D h_MC_Minv_pippim_All("h_MC_Minv_pippim_All", "", 200, 0., 3.);
    TH1D h_MC_Minv_pippim_RecSel("h_MC_Minv_pippim_RecSel", "", 200, 0., 3.);
    TH1D h_MC_Minv_pippim_Mx2Cut("h_MC_Minv_pippim_Mx2Cut", "", 200, 0., 3.);
    TH1D h_Minv_pippim_Mx2Cut_HasFT("h_Minv_pippim_Mx2Cut_HasFT", "", 200, 0., 3.);
    TH1D h_Minv_pippim_Mx2_Q2_cuts("h_Minv_pippim_Mx2_Q2_cuts", "", 200, 0., 3.);
    TH1D h_MMiss2_pippimprot("h_MMiss2_pippimprot", "", 200, -2., 2.);
    TH1D h_Q2_1("h_Q2_1", "", 400, 0., 1.2);
    TH1D h_Q2_HasFT1("h_Q2_HasFT1", "", 400, 0., 1.2);
    TH1D h_Q2_Mx2Cut("h_Q2_Mx2Cut", "", 400, 0., 1.2);
    TH1D h_Q2_Mx2Cut_HasFT("h_Q2_Mx2Cut_HasFT", "", 400, 0., 1.2);

    TH2D h_Q2_Mx2_1("h_Q2_Mx2_1", "", 200, -2., 2., 200, 0., 1.2);
    TH2D h_Q2_Mx2_With_emFT("h_Q2_Mx2_With_emFT", "", 200, -2., 2., 200, 0., 1.2);

    TH2D h_Minv_tM("h_Minv_tM", "", 200, 0., 3., 200, 0., 1.2);
    TH2D h_Minv_tM_MX2Cut("h_Minv_tM_MX2Cut", "", 200, 0., 3., 200, 0., 1.2);

    TH2D h_Q2_Eg1("h_Q2_Eg1", "", 200, 0., Eb, 200, 0., 0.2);
    TH2D h_Mx2_Eg1("h_Mx2_Eg1", "", 200, 0., Eb, 200, -2., 2.);

    TH1D h_n_emFT_0("h_n_emFT_0", "", 11, -0.5, 10.5);
    TH2D h_n_pip_pim_0("h_n_pip_pim_0", "", 11, -0.5, 10.5, 11, -0.5, 10.5);

    ROOT::Math::PxPyPzEVector L_beam(0., 0., Eb, Eb);
    ROOT::Math::PxPyPzEVector L_Targ(0., 0., 0., Mprot);

    hipo::reader reader;
    reader.open(inputFile.c_str());

    hipo::dictionary dict;
    reader.readDictionary(dict);

    RecEvent recEvent(dict);
    hipo::event event;

    MCEvent mcEvent(dict);

    long evCounter = 0;
    long nKept = 0;

    try {
        while (reader.next() == true) {
            reader.read(event);
            recEvent.Load(event);
            mcEvent.Load(event);
            evCounter = evCounter + 1;

            double MC_M_pipi = 0.;
            const bool hasMC_pipi = computeMC_Minv_pippim(mcEvent, MC_M_pipi);
            if (hasMC_pipi) {
                h_MC_Minv_pippim_All.Fill(MC_M_pipi);
            }

            std::vector<const RecPart*> v_pip, v_pim, v_prot, v_emFT;

            for (auto& p : recEvent.Particles()) {

                if (p.pid() == 211 && TMath::Abs(p.status()) >= STATUS_MIN && TMath::Abs(p.status()) < STATUS_MAX) {

                    h_vz_pip.Fill(p.vz());
                    h_chi2PID_pip.Fill(p.chi2pid());

                    double vt = 0.;
                    bool isVt = false;
                    if (computeVt(p, Mpi, vt)) {
                        h_vt_pip.Fill(vt);
                        isVt = vt > pip_vt_min && vt < pip_vt_max;
                    }

                    bool isVz = p.vz() > pip_vz_min && p.vz() < pip_vz_max;


                    if (isVt && isVz) {
                        v_pip.push_back(&p);
                    }

                } else if (p.pid() == -211 && TMath::Abs(p.status()) >= STATUS_MIN && TMath::Abs(p.status()) < STATUS_MAX) {

                    h_vz_pim.Fill(p.vz());
                    h_chi2PID_pim.Fill(p.chi2pid());

                    double vt = 0.;
                    bool isVt = false;
                    if (computeVt(p, Mpi, vt)) {
                        h_vt_pim.Fill(vt);
                        isVt = vt > pim_vt_min && vt < pim_vt_max;
                    }
                    bool isVz = p.vz() > pim_vz_min && p.vz() < pim_vz_max;

                    if (isVt && isVz) {
                        v_pim.push_back(&p);
                    }

                } else if (p.pid() == 2212) {

                    h_vz_prot.Fill(p.vz());
                    h_chi2PID_prot.Fill(p.chi2pid());

                    if (p.chi2pid() > protChi2PIDMin && p.chi2pid() < protChi2PIDMax) {

                        double vt = 0.;

                        if (computeVt(p, Mprot, vt)) {
                            h_vt_prot.Fill(vt);
                        } else {
                            cout << "No any scinitillator hit..." << endl;
                        }

                        bool isVz = p.vz() > prot_vz_Min && p.vz() < prot_vz_Max;
                        if ( isVz) {
                            v_prot.push_back(&p);
                        }
                    }
                }else if (p.pid() == 11 && TMath::Abs(p.status()) < 2000 ) {

                    v_emFT.push_back(&p);
                }
            }

            h_n_emFT_0.Fill(v_emFT.size());
            h_n_pip_pim_0.Fill(v_pim.size(), v_pip.size());

            if (v_pip.size() == 1 && v_pim.size() == 1 && v_prot.size() == 1) {
                nKept = nKept + 1;

                ROOT::Math::PxPyPzEVector L_pip, L_pim, L_prot, L_mis, L_gamma;
                L_pip.SetPxPyPzE(v_pip.at(0)->px(), v_pip.at(0)->py(), v_pip.at(0)->pz(), sqrt(v_pip.at(0)->p() * v_pip.at(0)->p() + Mpi * Mpi));
                L_pim.SetPxPyPzE(v_pim.at(0)->px(), v_pim.at(0)->py(), v_pim.at(0)->pz(), sqrt(v_pim.at(0)->p() * v_pim.at(0)->p() + Mpi * Mpi));
                L_prot.SetPxPyPzE(v_prot.at(0)->px(), v_prot.at(0)->py(), v_prot.at(0)->pz(), sqrt(v_prot.at(0)->p() * v_prot.at(0)->p() + Mprot * Mprot));

                L_mis = L_beam + L_Targ - L_pip - L_pim - L_prot;
                L_gamma = L_beam - L_mis;

                h_chi2PID_pip_FS.Fill(v_pip.at(0)->chi2pid());
                h_chi2PID_pim_FS.Fill(v_pim.at(0)->chi2pid());
                h_chi2PID_prot_FS.Fill(v_prot.at(0)->chi2pid());

                if (hasMC_pipi) {
                    h_MC_Minv_pippim_RecSel.Fill(MC_M_pipi);
                }

                double vt_FT;
                bool hasFT = false;
                const RecPart* p_emFT;
                if ( !v_emFT.empty() ) {
                    for (auto curPart : v_emFT) {
                        computeVtFT(*curPart, v_prot.at(0)->vz(), vt_FT);

                        vt_FT = vt_FT - v_prot.at(0)->vt();

                        h_vtem_FT.Fill(vt_FT);

                        if ( vt_FT > emFT_vt_min && vt_FT < emFT_vt_max ) {
                            hasFT = true;
                            p_emFT = curPart;
                        }
                    }

                }

                double M_pipi = (L_pip + L_pim).M();
                double Mx2 = L_mis.M2();
                double Q2 = 2*Eb*L_mis.P()*(1 - cos(L_mis.Theta()));
                double Eg = L_gamma.E();
                double tM = 2*Mprot*(L_prot.E() - Mprot);

                h_Minv_tM.Fill(M_pipi, tM);
                h_Minv_pippim.Fill(M_pipi);
                h_MMiss2_pippimprot.Fill(Mx2);
                h_Q2_1.Fill(Q2);
                h_Q2_Mx2_1.Fill(Mx2, Q2);

                h_Q2_Eg1.Fill(Eg, Q2);
                h_Mx2_Eg1.Fill(Eg, Mx2);

                if (hasFT) {
                    h_Q2_Mx2_With_emFT.Fill(Mx2, Q2);
                    h_Q2_HasFT1.Fill(Q2);
                }

                if ( Mx2 > Mx2_Min && Mx2 < Mx2_Max ) {
                    h_Minv_pippim_Mx2Cut.Fill(M_pipi);

                    if (hasMC_pipi) {
                        h_MC_Minv_pippim_Mx2Cut.Fill(MC_M_pipi);
                    }

                    h_Q2_Mx2Cut.Fill(Q2);
                    h_Minv_tM_MX2Cut.Fill(M_pipi, tM);

                    if (hasFT) {
                        h_Minv_pippim_Mx2Cut_HasFT.Fill(M_pipi);
                        h_Q2_Mx2Cut_HasFT.Fill(Q2);
                    }

                    if ( Q2 < Q2_Max ) {
                        h_Minv_pippim_Mx2_Q2_cuts.Fill(M_pipi);
                    }
                }

            }

            if (evCounter % 50000 == 0) {
                cout.flush() << "Processed " << evCounter << " events, kept " << nKept << " \r";
            }
        }
    } catch (std::exception& e) {
        cerr << e.what() << endl;
    }

    cout << endl << "Done. Processed " << evCounter << " events, kept " << nKept
            << " -> " << outputFile << endl;

    gDirectory->Write();
    file_out.Close();

    return 0;
}
