/*
 * File:   CompareHists.cc
 * Author: rafopar
 *
 * Compares every histogram present under the same name in two root files,
 * overlaying 1D histograms (each normalized to a maximum of 1) and, for 2D
 * histograms, laying out a 2x2 canvas with the raw 2D histogram from each
 * file on the left and normalized 1D X/Y projections overlaid on the right.
 * Every comparison is written as one page of a single multi-page pdf.
 *
 * Usage:
 *   CompareHists.exe <file1.root> <file2.root> <keyword1> <keyword2>
 *
 * <keyword1>/<keyword2> label the legends (file1 in red, file2 in blue) and
 * name the output file:
 *   ./Figs/RhoTail_Comparisons_<keyword1>_<keyword2>.pdf
 *
 * Titles and axis labels are not read from the histograms -- AnaPiPiProt.cc
 * creates every histogram with an empty title, so this file keeps its own
 * name -> "Title;XTitle;YTitle" table (HistTitles) instead, matched against
 * the selection logic in RhoTail/AnaPiPiProt.cc. Cut values applied there
 * (VCuts/HCuts) are drawn as dashed vertical/horizontal reference lines.
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

using namespace std;

namespace {

    // Scales h in place so that its maximum bin content becomes 1.
    void NormalizeToUnitMax(TH1* h) {
        const double max = h->GetMaximum();
        if (max > 0.) {
            h->Scale(1. / max);
        }
    }

    void StyleAsFile1(TH1* h) {
        h->SetLineColor(kRed);
        h->SetLineWidth(2);
        h->SetMarkerColor(kRed);
    }

    void StyleAsFile2(TH1* h) {
        h->SetLineColor(kBlue);
        h->SetLineWidth(2);
        h->SetMarkerColor(kBlue);
    }

    // Histogram titles, as "Title;XAxisTitle;YAxisTitle" (ROOT's TH1::SetTitle
    // parses this the same way the constructor does). Keyed by name against
    // what RhoTail/AnaPiPiProt.cc fills each histogram with.
    const unordered_map<string, string> HistTitles = {
        {"h_vz_pip", "V_{z} of #pi^{+};V_{z}^{#pi^{+}} (cm);Counts"},
        {"h_vz_pim", "V_{z} of #pi^{-};V_{z}^{#pi^{-}} (cm);Counts"},
        {"h_vz_prot", "V_{z} of proton;V_{z}^{p} (cm);Counts"},
        {"h_chi2PID_pip", "#chi^{2}_{PID} of #pi^{+};#chi^{2}_{PID}(#pi^{+});Counts"},
        {"h_chi2PID_pim", "#chi^{2}_{PID} of #pi^{-};#chi^{2}_{PID}(#pi^{-});Counts"},
        {"h_chi2PID_prot", "#chi^{2}_{PID} of proton;#chi^{2}_{PID}(p);Counts"},
        {"h_vt_pip", "Vertex time, #pi^{+};#Deltat_{v}^{#pi^{+}} (ns);Counts"},
        {"h_vt_pim", "Vertex time, #pi^{-};#Deltat_{v}^{#pi^{-}} (ns);Counts"},
        {"h_vt_prot", "Vertex time, proton;#Deltat_{v}^{p} (ns);Counts"},
        {"h_vtem_FT", "Vertex time, FT electron;#Deltat_{v}^{e,FT} (ns);Counts"},
        {"h_Minv_pippim", "M(#pi^{+}#pi^{-}), no cuts;M(#pi^{+}#pi^{-}) (GeV);Counts"},
        {"h_Minv_pippim_Mx2Cut", "M(#pi^{+}#pi^{-}), M_{x}^{2} cut;M(#pi^{+}#pi^{-}) (GeV);Counts"},
        {"h_Minv_pippim_Mx2Cut_HasFT", "M(#pi^{+}#pi^{-}), M_{x}^{2} cut, FT e^{-};M(#pi^{+}#pi^{-}) (GeV);Counts"},
        {"h_Minv_pippim_Mx2_Q2_cuts", "M(#pi^{+}#pi^{-}), M_{x}^{2} & Q^{2} cuts;M(#pi^{+}#pi^{-}) (GeV);Counts"},
        {"h_MMiss2_pippimprot", "Missing mass squared, ep#rightarrowe#pi^{+}#pi^{-}pX;M_{x}^{2} (GeV^{2});Counts"},
        {"h_Q2_1", "Q^{2}, no cuts;Q^{2} (GeV^{2});Counts"},
        {"h_Q2_HasFT1", "Q^{2}, FT electron tagged;Q^{2} (GeV^{2});Counts"},
        {"h_Q2_Mx2Cut", "Q^{2}, M_{x}^{2} cut;Q^{2} (GeV^{2});Counts"},
        {"h_Q2_Mx2Cut_HasFT", "Q^{2}, M_{x}^{2} cut, FT e^{-};Q^{2} (GeV^{2});Counts"},
        {"h_Q2_Mx2_1", "Q^{2} vs M_{x}^{2};M_{x}^{2} (GeV^{2});Q^{2} (GeV^{2})"},
        {"h_Q2_Mx2_With_emFT", "Q^{2} vs M_{x}^{2}, FT e^{-} tagged;M_{x}^{2} (GeV^{2});Q^{2} (GeV^{2})"},
        {"h_Minv_tM", "M(#pi^{+}#pi^{-}) vs t;M(#pi^{+}#pi^{-}) (GeV);t (GeV^{2})"},
        {"h_Minv_tM_MX2Cut", "M(#pi^{+}#pi^{-}) vs t, M_{x}^{2} cut;M(#pi^{+}#pi^{-}) (GeV);t (GeV^{2})"},
        {"h_Q2_Eg1", "Q^{2} vs E_{#gamma};E_{#gamma} (GeV);Q^{2} (GeV^{2})"},
        {"h_Mx2_Eg1", "M_{x}^{2} vs E_{#gamma};E_{#gamma} (GeV);M_{x}^{2} (GeV^{2})"},
        {"h_n_emFT_0", "Number of FT electron candidates;N_{e}^{FT};Counts"},
    };

    // Cut values from RhoTail/AnaPiPiProt.cc, keyed by histogram name and
    // drawn as dashed reference lines: VCuts on the histogram's X axis (every
    // 1D histogram, plus the X variable of a 2D one), HCuts on a 2D
    // histogram's Y axis only.
    const unordered_map<string, vector<double> > VCuts = {
        {"h_vz_pip", {-10., 2.}}, // pip_vz_min, pip_vz_max
        {"h_vz_pim", {-8., 3.}}, // pim_vz_min, pim_vz_max
        {"h_vz_prot", {-7., 1.}}, // prot_vz_Min, prot_vz_Max
        {"h_chi2PID_prot", {-10., 10.}}, // protChi2PIDMin, protChi2PIDMax
        {"h_vt_pip", {-0.8, 0.8}}, // pip_vt_min, pip_vt_max
        {"h_vt_pim", {-0.8, 0.8}}, // pim_vt_min, pim_vt_max
        {"h_vtem_FT", {-2., 2.}}, // emFT_vt_min, emFT_vt_max
        {"h_MMiss2_pippimprot", {-0.5, 0.5}}, // Mx2_Min, Mx2_Max
        {"h_Q2_1", {0.15}}, // Q2_Max
        {"h_Q2_HasFT1", {0.15}}, // Q2_Max
        {"h_Q2_Mx2Cut", {0.15}}, // Q2_Max
        {"h_Q2_Mx2Cut_HasFT", {0.15}}, // Q2_Max
        {"h_Q2_Mx2_1", {-0.5, 0.5}}, // Mx2_Min, Mx2_Max (X axis)
        {"h_Q2_Mx2_With_emFT", {-0.5, 0.5}}, // Mx2_Min, Mx2_Max (X axis)
    };

    const unordered_map<string, vector<double> > HCuts = {
        {"h_Q2_Mx2_1", {0.15}}, // Q2_Max (Y axis)
        {"h_Q2_Mx2_With_emFT", {0.15}}, // Q2_Max (Y axis)
        {"h_Q2_Eg1", {0.15}}, // Q2_Max (Y axis)
        {"h_Mx2_Eg1", {-0.5, 0.5}}, // Mx2_Min, Mx2_Max (Y axis)
    };

    const vector<double> kNoCuts;

    string TitleFor(const string& name) {
        auto it = HistTitles.find(name);
        return (it != HistTitles.end()) ? it->second : name;
    }

    const vector<double>& VCutsFor(const string& name) {
        auto it = VCuts.find(name);
        return (it != VCuts.end()) ? it->second : kNoCuts;
    }

    const vector<double>& HCutsFor(const string& name) {
        auto it = HCuts.find(name);
        return (it != HCuts.end()) ? it->second : kNoCuts;
    }

    // Splits a "Title;XTitle;YTitle" string into its (up to) 3 parts, so a 2D
    // histogram's projections can carry over the matching axis title.
    vector<string> SplitTitle(const string& t) {
        vector<string> parts;
        size_t start = 0;
        for (int i = 0; i < 2 && start != string::npos; ++i) {
            size_t pos = t.find(';', start);
            if (pos == string::npos) {
                parts.push_back(t.substr(start));
                start = string::npos;
            } else {
                parts.push_back(t.substr(start, pos - start));
                start = pos + 1;
            }
        }
        if (start != string::npos) {
            parts.push_back(t.substr(start));
        }
        while (parts.size() < 3) {
            parts.push_back("");
        }
        return parts;
    }

    // Inserts " (suffix)" into the title part of a "Title;X;Y" string, leaving
    // the axis titles themselves untouched.
    string WithSuffix(const string& title, const string& suffix) {
        size_t pos = title.find(';');
        const string mainTitle = (pos == string::npos) ? title : title.substr(0, pos);
        const string rest = (pos == string::npos) ? string() : title.substr(pos);
        return mainTitle + " (" + suffix + ")" + rest;
    }

    // Dashed vertical lines at each cut value in `xs`, spanning h's current Y
    // range (as set by DrawOverlay via SetMinimum/SetMaximum).
    void DrawVerticalCuts(TH1* h, const vector<double>& xs) {
        for (double x : xs) {
            TLine* l = new TLine(x, h->GetMinimum(), x, h->GetMaximum());
            l->SetLineColor(kBlack);
            l->SetLineStyle(2);
            l->Draw();
        }
    }

    // Dashed vertical/horizontal lines spanning h's full axis range, for cuts
    // on a 2D histogram's X and Y variables respectively.
    void DrawCuts2D(TH2* h, const vector<double>& xs, const vector<double>& ys) {
        for (double x : xs) {
            TLine* l = new TLine(x, h->GetYaxis()->GetXmin(), x, h->GetYaxis()->GetXmax());
            l->SetLineColor(kBlack);
            l->SetLineStyle(2);
            l->Draw();
        }
        for (double y : ys) {
            TLine* l = new TLine(h->GetXaxis()->GetXmin(), y, h->GetXaxis()->GetXmax(), y);
            l->SetLineColor(kBlack);
            l->SetLineStyle(2);
            l->Draw();
        }
    }

    // Names of every TH1-derived object at the top level of f, in key order,
    // de-duplicated (a file can list multiple key cycles for the same name).
    vector<string> HistNames(TFile* f) {
        vector<string> names;
        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*) next())) {
            TClass* cl = TClass::GetClass(key->GetClassName());
            if (cl == nullptr || !cl->InheritsFrom("TH1")) {
                continue;
            }
            if (find(names.begin(), names.end(), string(key->GetName())) == names.end()) {
                names.push_back(key->GetName());
            }
        }
        return names;
    }

    // Draws h1 and h2 (already styled/normalized) on the current pad,
    // overlaid, with a legend labeling them by keyword.
    void DrawOverlay(TH1* h1, TH1* h2, const string& keyword1, const string& keyword2,
            const string& title, const vector<double>& cutLines = {}) {
        const double ymax = 1.15 * max(h1->GetMaximum(), h2->GetMaximum());
        h1->SetTitle(title.c_str());
        h1->SetMinimum(0.);
        h1->SetMaximum(ymax > 0. ? ymax : 1.);
        h1->Draw("HIST");
        h2->Draw("HIST SAME");

        DrawVerticalCuts(h1, cutLines);

        TLegend* leg = new TLegend(0.65, 0.78, 0.89, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h1, keyword1.c_str(), "l");
        leg->AddEntry(h2, keyword2.c_str(), "l");
        leg->Draw();
    }
}

int main(int argc, char** argv) {

    if (argc != 5) {
        cerr << "Usage: " << argv[0]
                << " <file1.root> <file2.root> <keyword1> <keyword2>" << endl;
        return 1;
    }

    const string file1Name = argv[1];
    const string file2Name = argv[2];
    const string keyword1 = argv[3];
    const string keyword2 = argv[4];

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TFile* f1 = TFile::Open(file1Name.c_str(), "READ");
    if (f1 == nullptr || f1->IsZombie()) {
        cerr << "Could not open " << file1Name << endl;
        return 1;
    }
    TFile* f2 = TFile::Open(file2Name.c_str(), "READ");
    if (f2 == nullptr || f2->IsZombie()) {
        cerr << "Could not open " << file2Name << endl;
        return 1;
    }

    // Match by name against file1's histogram list; warn about anything in
    // file1 that has no counterpart in file2.
    const vector<string> names1 = HistNames(f1);
    vector<string> commonNames;
    for (const string& name : names1) {
        if (f2->Get(name.c_str()) != nullptr) {
            commonNames.push_back(name);
        } else {
            cerr << "Warning: '" << name << "' found in " << file1Name
                    << " but not in " << file2Name << " -- skipping." << endl;
        }
    }

    if (commonNames.empty()) {
        cerr << "No histograms with matching names found in both files." << endl;
        return 1;
    }

    gSystem->mkdir("./Figs", true);
    const TString outPdf = Form("./Figs/RhoTail_Comparisons_%s_%s.pdf",
            keyword1.c_str(), keyword2.c_str());

    TCanvas c1D("c1D", "c1D", 900, 700);
    TCanvas c2D("c2D", "c2D", 1100, 900);

    for (size_t i = 0; i < commonNames.size(); ++i) {
        const string& name = commonNames[i];
        const string pageOpt = (i == 0) ? "(" : (i == commonNames.size() - 1) ? ")" : "";

        TH1* h1 = dynamic_cast<TH1*> (f1->Get(name.c_str()));
        TH1* h2 = dynamic_cast<TH1*> (f2->Get(name.c_str()));
        if (h1 == nullptr || h2 == nullptr) {
            continue;
        }

        if (!h1->InheritsFrom("TH2")) {
            // ---- Plain 1D comparison ----
            TH1* c1 = (TH1*) h1->Clone((name + "_cmp1").c_str());
            TH1* c2 = (TH1*) h2->Clone((name + "_cmp2").c_str());
            c1->SetDirectory(nullptr);
            c2->SetDirectory(nullptr);
            NormalizeToUnitMax(c1);
            NormalizeToUnitMax(c2);
            StyleAsFile1(c1);
            StyleAsFile2(c2);

            c1D.cd();
            c1D.Clear();
            DrawOverlay(c1, c2, keyword1, keyword2, TitleFor(name), VCutsFor(name));
            c1D.Print(outPdf + pageOpt, "pdf");

            delete c1;
            delete c2;
        } else {
            // ---- 2D: raw histograms on the left, normalized projections on the right ----
            TH2* h1_2d = (TH2*) h1->Clone((name + "_2d1").c_str());
            TH2* h2_2d = (TH2*) h2->Clone((name + "_2d2").c_str());
            h1_2d->SetDirectory(nullptr);
            h2_2d->SetDirectory(nullptr);

            const string title2D = TitleFor(name);
            const vector<string> titleParts = SplitTitle(title2D);
            const string& xTitle = titleParts[1];
            const string& yTitle = titleParts[2];
            const vector<double>& xCuts = VCutsFor(name);
            const vector<double>& yCuts = HCutsFor(name);

            c2D.Clear();
            c2D.Divide(2, 2);

            c2D.cd(1);
            gPad->SetRightMargin(0.15);
            h1_2d->SetTitle(WithSuffix(title2D, keyword1).c_str());
            h1_2d->Draw("COLZ");
            DrawCuts2D(h1_2d, xCuts, yCuts);

            c2D.cd(3);
            gPad->SetRightMargin(0.15);
            h2_2d->SetTitle(WithSuffix(title2D, keyword2).c_str());
            h2_2d->Draw("COLZ");
            DrawCuts2D(h2_2d, xCuts, yCuts);

            TH1D* px1 = h1_2d->ProjectionX((name + "_px1").c_str());
            TH1D* px2 = h2_2d->ProjectionX((name + "_px2").c_str());
            px1->SetDirectory(nullptr);
            px2->SetDirectory(nullptr);
            NormalizeToUnitMax(px1);
            NormalizeToUnitMax(px2);
            StyleAsFile1(px1);
            StyleAsFile2(px2);

            const string pxTitle = (xTitle.empty() ? (name + " : Projection X") : (xTitle + " projection"))
                    + ";" + xTitle + ";Normalized counts";

            c2D.cd(2);
            DrawOverlay(px1, px2, keyword1, keyword2, pxTitle, xCuts);

            TH1D* py1 = h1_2d->ProjectionY((name + "_py1").c_str());
            TH1D* py2 = h2_2d->ProjectionY((name + "_py2").c_str());
            py1->SetDirectory(nullptr);
            py2->SetDirectory(nullptr);
            NormalizeToUnitMax(py1);
            NormalizeToUnitMax(py2);
            StyleAsFile1(py1);
            StyleAsFile2(py2);

            const string pyTitle = (yTitle.empty() ? (name + " : Projection Y") : (yTitle + " projection"))
                    + ";" + yTitle + ";Normalized counts";

            c2D.cd(4);
            // Cuts on the 2D's Y variable become vertical lines here, since
            // this projection's X axis IS that Y variable.
            DrawOverlay(py1, py2, keyword1, keyword2, pyTitle, yCuts);

            c2D.Print(outPdf + pageOpt, "pdf");

            delete px1;
            delete px2;
            delete py1;
            delete py2;
            delete h1_2d;
            delete h2_2d;
        }
    }

    cout << "Wrote " << commonNames.size() << " comparison page(s) to "
            << outPdf << endl;

    f1->Close();
    f2->Close();

    return 0;
}
