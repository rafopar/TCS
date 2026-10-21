/*
 * File:   RhoFitGUI.cc
 * Author: rafopar
 *
 * Interactive ROOT GUI to fit the pi-pi+ invariant mass distribution with
 * the sum of three relativistic Breit-Wigners (rho, f2(1270), f0(1500)),
 * optionally convolved with a Gaussian to account for the mass resolution.
 *
 * Run interpreted or ACLiC-compiled from the RhoTail directory:
 *   root -l RhoFitGUI.cc
 *   root -l RhoFitGUI.cc+
 *
 * By default it opens Hists/Hists_PionPairSkim_All.root (as produced by
 * RunAnaPiPiProt.py's hadd merge); pass a different path as the argument if
 * needed:
 *   root -l 'RhoFitGUI.cc("some/other/file.root")'
 *
 * Use the combo box to pick which histogram to fit -- e.g.
 * h_Minv_pippim_Mx2_Q2_cuts.
 */

#include <cmath>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1D.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TApplication.h>

#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGComboBox.h>
#include <TGNumberEntry.h>
#include <TRootEmbeddedCanvas.h>

// ============================== Physics ====================================

namespace {
    const double Mpi = 0.13957; // GeV, charged pion mass

    // Two-body breakup momentum of a resonance of mass m decaying into two
    // pions of mass Mpi, in the resonance rest frame.
    double Q2Body(double m) {
        if (m < 2. * Mpi) {
            return 0.;
        }
        return std::sqrt(m * m / 4. - Mpi * Mpi);
    }

    // Relativistic Breit-Wigner with a mass-dependent width (Jackson form)
    // for a resonance of spin/orbital-angular-momentum L decaying into two
    // identical pions.
    double RelBW(double m, double M0, double Gamma0, double L) {
        double q = Q2Body(m);
        double q0 = Q2Body(M0);
        if (q0 <= 0. || m <= 0.) {
            return 0.;
        }

        double GammaM = Gamma0 * std::pow(q / q0, 2. * L + 1.) * (M0 / m);
        double num = M0 * GammaM;
        double den = (m * m - M0 * M0) * (m * m - M0 * M0) + M0 * M0 * GammaM * GammaM;

        return num / den;
    }

    // Sum of three relativistic Breit-Wigners: rho (P-wave, L=1), f2(1270)
    // (D-wave, L=2), and f0(1500) (S-wave, L=0). Parameters: [0]=A_rho,
    // [1]=M_rho, [2]=Gamma_rho, [3]=A_f2, [4]=M_f2, [5]=Gamma_f2, [6]=A_f0,
    // [7]=M_f0, [8]=Gamma_f0.
    double SumRBW(double* x, double* par) {
        double m = x[0];

        double bwRho = par[0] * RelBW(m, par[1], par[2], 1.);
        double bwF2 = par[3] * RelBW(m, par[4], par[5], 2.);
        double bwF0 = par[6] * RelBW(m, par[7], par[8], 0.);

        return bwRho + bwF2 + bwF0;
    }

    // The individual rho (P-wave), f2 (D-wave), and f0 (S-wave) components,
    // for drawing each resonance separately. Parameters: [0]=A, [1]=M,
    // [2]=Gamma.
    double RhoOnly(double* x, double* par) {
        return par[0] * RelBW(x[0], par[1], par[2], 1.);
    }

    double F2Only(double* x, double* par) {
        return par[0] * RelBW(x[0], par[1], par[2], 2.);
    }

    double F0Only(double* x, double* par) {
        return par[0] * RelBW(x[0], par[1], par[2], 0.);
    }

    const int kNpxSmooth = 4500;

    // Normalized Gaussian resolution kernel, mean fixed at zero, sigma is
    // the single free parameter [0].
    double GausKernel(double* x, double* par) {
        double sigma = par[0];
        if (sigma <= 0.) {
            return (x[0] == 0.) ? 1. : 0.;
        }
        return 1. / (sigma * std::sqrt(2. * TMath::Pi())) * std::exp(-0.5 * (x[0] / sigma) * (x[0] / sigma));
    }

    const int kNParBase = 9;
    const char* kParNames[kNParBase + 1] = {
        "A(rho)", "M(rho)", "#Gamma(rho)",
        "A(f2)", "M(f2)", "#Gamma(f2)",
        "A(f0)", "M(f0)", "#Gamma(f0)",
        "#sigma(res)"
    };
    const double kParInit[kNParBase + 1] = {
        1000., 0.770, 0.150,
        100., 1.270, 0.185,
        50., 1.506, 0.112,
        0.020
    };
    const double kParMin[kNParBase + 1] = {
        0., 0.6, 0.05,
        0., 1.1, 0.05,
        0., 1.3, 0.05,
        0.001
    };
    const double kParMax[kNParBase + 1] = {
        1.e9, 0.9, 0.3,
        1.e9, 1.4, 0.4,
        1.e9, 1.7, 0.3,
        0.1
    };
}

// =============================== GUI class ==================================

struct ParamRow {
    TGLabel* label = nullptr;
    TGNumberEntry* value = nullptr;
    TGCheckButton* fix = nullptr;
    TGNumberEntry* min = nullptr;
    TGNumberEntry* max = nullptr;
};

class RhoFitPanel : public TGMainFrame {
public:
    RhoFitPanel(const TGWindow* p, UInt_t w, UInt_t h, const char* fileName);
    ~RhoFitPanel() override;

    // Slots
    void OnHistSelected(Int_t id);
    void OnConvolveToggled();
    void DoFit();
    void DoDrawFunction();
    void DoSave();
    void DoExit();

private:
    TFile* fFile = nullptr;
    TH1D* fHist = nullptr;
    TF1* fFitFunc = nullptr; // total (rho + f2 + f0) function used for fitting/drawing
    TF1* fRhoFunc = nullptr; // rho-only component, drawn separately
    TF1* fF2Func = nullptr;  // f2-only component, drawn separately
    TF1* fF0Func = nullptr;  // f0-only component, drawn separately
    TLegend* fLegend = nullptr;

    TRootEmbeddedCanvas* fEmbCanvas = nullptr;
    TGComboBox* fHistCombo = nullptr;

    TGNumberEntry* fFitMinEntry = nullptr;
    TGNumberEntry* fFitMaxEntry = nullptr;
    TGNumberEntry* fDrawMinEntry = nullptr;
    TGNumberEntry* fDrawMaxEntry = nullptr;

    TGCheckButton* fConvolveCheck = nullptr;
    TGLabel* fStatusLabel = nullptr;

    std::vector<ParamRow> fParamRows; // 10 rows: 9 base (3 resonances x 3) + 1 resolution sigma
    std::vector<std::string> fHistNames;

    void BuildGUI();
    void PopulateHistList();
    void LoadHistogram(const char* name);
    TF1* BuildFunction(bool convolved, double xmin, double xmax) const;
    TF1* BuildComponentFunction(int which, bool convolved, double xmin, double xmax) const;
    void ApplyWidgetsToFunction(TF1* f) const;
    void RebuildComponents(bool convolved, double xmin, double xmax);
    void SyncComponentParams();
    void Redraw();
};

RhoFitPanel::RhoFitPanel(const TGWindow* p, UInt_t w, UInt_t h, const char* fileName)
: TGMainFrame(p, w, h) {

    fFile = TFile::Open(fileName, "READ");
    if (!fFile || fFile->IsZombie()) {
        Error("RhoFitPanel", "Could not open file \"%s\"", fileName);
    }

    BuildGUI();
    PopulateHistList();

    SetWindowName("Rho/f2/f0 Invariant Mass Fitter");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}

RhoFitPanel::~RhoFitPanel() {
    delete fRhoFunc;
    delete fF2Func;
    delete fF0Func;
    delete fFitFunc;
    delete fLegend;
    if (fFile) {
        fFile->Close();
    }
}

void RhoFitPanel::BuildGUI() {

    auto* mainHorz = new TGHorizontalFrame(this, 1, 1);

    // ----- Left: canvas -----
    fEmbCanvas = new TRootEmbeddedCanvas("fEmbCanvas", mainHorz, 800, 650);
    mainHorz->AddFrame(fEmbCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY | kLHintsLeft, 5, 5, 5, 5));

    // ----- Right: controls -----
    auto* controls = new TGVerticalFrame(mainHorz, 380, 650);

    // Histogram selector
    auto* histFrame = new TGGroupFrame(controls, "Histogram", kVerticalFrame);
    fHistCombo = new TGComboBox(histFrame, -1);
    fHistCombo->Resize(340, 24);
    fHistCombo->Connect("Selected(Int_t)", "RhoFitPanel", this, "OnHistSelected(Int_t)");
    histFrame->AddFrame(fHistCombo, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));
    controls->AddFrame(histFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));

    // Fit / draw ranges
    auto* rangeFrame = new TGGroupFrame(controls, "Ranges (GeV)", kVerticalFrame);

    auto* fitRangeRow = new TGHorizontalFrame(rangeFrame, 1, 1);
    fitRangeRow->AddFrame(new TGLabel(fitRangeRow, "Fit range:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 6, 2, 2));
    fFitMinEntry = new TGNumberEntry(fitRangeRow, 0.4, 7, -1, TGNumberFormat::kNESReal);
    fitRangeRow->AddFrame(fFitMinEntry, new TGLayoutHints(kLHintsLeft, 2, 6, 2, 2));
    fFitMaxEntry = new TGNumberEntry(fitRangeRow, 1.6, 7, -1, TGNumberFormat::kNESReal);
    fitRangeRow->AddFrame(fFitMaxEntry, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    rangeFrame->AddFrame(fitRangeRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    auto* drawRangeRow = new TGHorizontalFrame(rangeFrame, 1, 1);
    drawRangeRow->AddFrame(new TGLabel(drawRangeRow, "Draw range:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 2, 2));
    fDrawMinEntry = new TGNumberEntry(drawRangeRow, 0.3, 7, -1, TGNumberFormat::kNESReal);
    drawRangeRow->AddFrame(fDrawMinEntry, new TGLayoutHints(kLHintsLeft, 2, 6, 2, 2));
    fDrawMaxEntry = new TGNumberEntry(drawRangeRow, 1.8, 7, -1, TGNumberFormat::kNESReal);
    drawRangeRow->AddFrame(fDrawMaxEntry, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    rangeFrame->AddFrame(drawRangeRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    controls->AddFrame(rangeFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));

    // Parameters table
    auto* parFrame = new TGGroupFrame(controls, "Parameters (fix / range)", kVerticalFrame);

    auto* headerRow = new TGHorizontalFrame(parFrame, 1, 1);
    headerRow->AddFrame(new TGLabel(headerRow, "                "), new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    headerRow->AddFrame(new TGLabel(headerRow, "  value  "), new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    headerRow->AddFrame(new TGLabel(headerRow, "fix"), new TGLayoutHints(kLHintsLeft, 8, 2, 2, 2));
    headerRow->AddFrame(new TGLabel(headerRow, "   min   "), new TGLayoutHints(kLHintsLeft, 8, 2, 2, 2));
    headerRow->AddFrame(new TGLabel(headerRow, "   max   "), new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    parFrame->AddFrame(headerRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    fParamRows.resize(kNParBase + 1);
    for (int i = 0; i < kNParBase + 1; i++) {
        auto* row = new TGHorizontalFrame(parFrame, 1, 1);
        ParamRow& pr = fParamRows[i];

        pr.label = new TGLabel(row, kParNames[i]);
        row->AddFrame(pr.label, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 2, 2));

        pr.value = new TGNumberEntry(row, kParInit[i], 7, -1, TGNumberFormat::kNESReal);
        pr.value->Resize(70, pr.value->GetDefaultHeight());
        row->AddFrame(pr.value, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

        pr.fix = new TGCheckButton(row, "");
        row->AddFrame(pr.fix, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 2, 2));

        pr.min = new TGNumberEntry(row, kParMin[i], 7, -1, TGNumberFormat::kNESReal);
        pr.min->Resize(60, pr.min->GetDefaultHeight());
        row->AddFrame(pr.min, new TGLayoutHints(kLHintsLeft, 8, 2, 2, 2));

        pr.max = new TGNumberEntry(row, kParMax[i], 7, -1, TGNumberFormat::kNESReal);
        pr.max->Resize(60, pr.max->GetDefaultHeight());
        row->AddFrame(pr.max, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

        parFrame->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    }
    controls->AddFrame(parFrame, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));

    // Convolution toggle
    fConvolveCheck = new TGCheckButton(controls, "Convolve with Gaussian (mass resolution)");
    fConvolveCheck->Connect("Clicked()", "RhoFitPanel", this, "OnConvolveToggled()");
    controls->AddFrame(fConvolveCheck, new TGLayoutHints(kLHintsExpandX, 8, 4, 8, 4));
    // sigma(res) row (last row of fParamRows) is only meaningful with convolution on
    fParamRows[kNParBase].value->SetState(kFALSE);
    fParamRows[kNParBase].min->SetState(kFALSE);
    fParamRows[kNParBase].max->SetState(kFALSE);
    fParamRows[kNParBase].fix->SetState(kButtonDisabled);

    // Action buttons
    auto* buttonRow = new TGHorizontalFrame(controls, 1, 1);
    auto* fitButton = new TGTextButton(buttonRow, "Fit");
    fitButton->Connect("Clicked()", "RhoFitPanel", this, "DoFit()");
    buttonRow->AddFrame(fitButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 4, 2, 2));

    auto* drawButton = new TGTextButton(buttonRow, "Draw Function");
    drawButton->Connect("Clicked()", "RhoFitPanel", this, "DoDrawFunction()");
    buttonRow->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 4, 2, 2));

    auto* saveButton = new TGTextButton(buttonRow, "Save");
    saveButton->Connect("Clicked()", "RhoFitPanel", this, "DoSave()");
    buttonRow->AddFrame(saveButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 2, 2, 2));

    controls->AddFrame(buttonRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 4));

    // Exit on its own row, separated from the routine actions above.
    auto* exitButton = new TGTextButton(controls, "Exit");
    exitButton->Connect("Clicked()", "RhoFitPanel", this, "DoExit()");
    controls->AddFrame(exitButton, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    // Status line
    fStatusLabel = new TGLabel(controls, "Ready.");
    controls->AddFrame(fStatusLabel, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 8));

    mainHorz->AddFrame(controls, new TGLayoutHints(kLHintsExpandY | kLHintsLeft, 4, 4, 4, 4));

    AddFrame(mainHorz, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
}

void RhoFitPanel::PopulateHistList() {
    fHistCombo->RemoveAll();
    fHistNames.clear();

    if (!fFile) {
        return;
    }

    int id = 0;
    int defaultId = 0;
    TIter next(fFile->GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next())) {
        TClass* cl = TClass::GetClass(key->GetClassName());
        if (!cl || !cl->InheritsFrom("TH1")) {
            continue;
        }

        const char* name = key->GetName();
        fHistCombo->AddEntry(name, id);
        fHistNames.push_back(name);

        if (std::string(name) == "h_Minv_pippim_Mx2_Q2_cuts") {
            defaultId = id;
        }

        id++;
    }

    if (id > 0) {
        fHistCombo->Select(defaultId, kTRUE);
    }
}

void RhoFitPanel::LoadHistogram(const char* name) {
    if (!fFile) {
        return;
    }

    TH1D* h = dynamic_cast<TH1D*>(fFile->Get(name));
    if (!h) {
        Error("LoadHistogram", "Histogram \"%s\" not found", name);
        return;
    }

    fHist = h;
    fHist->SetDirectory(nullptr);

    // Drop any function drawn for the previously selected histogram.
    delete fFitFunc; fFitFunc = nullptr;
    delete fRhoFunc; fRhoFunc = nullptr;
    delete fF2Func; fF2Func = nullptr;
    delete fF0Func; fF0Func = nullptr;

    fFitMinEntry->SetNumber(fHist->GetXaxis()->GetXmin());
    fFitMaxEntry->SetNumber(fHist->GetXaxis()->GetXmax());
    fDrawMinEntry->SetNumber(fHist->GetXaxis()->GetXmin());
    fDrawMaxEntry->SetNumber(fHist->GetXaxis()->GetXmax());

    Redraw();
}

TF1* RhoFitPanel::BuildFunction(bool convolved, double xmin, double xmax) const {
    if (!convolved) {
        auto* f = new TF1("fSumRBW", SumRBW, xmin, xmax, kNParBase);
        for (int i = 0; i < kNParBase; i++) {
            f->SetParName(i, kParNames[i]);
        }
        f->SetNpx(kNpxSmooth);
        return f;
    }

    auto* fBase = new TF1("fBaseRBW", SumRBW, xmin, xmax, kNParBase);
    auto* fKernel = new TF1("fGausKernel", GausKernel, -1., 1., 1);

    TF1Convolution conv(fBase, fKernel, xmin, xmax, true);
    conv.SetNofPointsFFT(1000);

    auto* f = new TF1("fConvRBW", conv, xmin, xmax, conv.GetNpar());
    for (int i = 0; i < kNParBase + 1; i++) {
        f->SetParName(i, kParNames[i]);
    }
    f->SetNpx(kNpxSmooth);
    return f;
}

TF1* RhoFitPanel::BuildComponentFunction(int which, bool convolved, double xmin, double xmax) const {
    double (*func)(double*, double*) = (which == 0) ? RhoOnly : (which == 1) ? F2Only : F0Only;
    const char* baseName = (which == 0) ? "fRhoBase" : (which == 1) ? "fF2Base" : "fF0Base";
    const char* kernelName = (which == 0) ? "fRhoKernel" : (which == 1) ? "fF2Kernel" : "fF0Kernel";
    const char* fullName = (which == 0) ? "fRhoComp" : (which == 1) ? "fF2Comp" : "fF0Comp";

    if (!convolved) {
        auto* f = new TF1(fullName, func, xmin, xmax, 3);
        f->SetNpx(kNpxSmooth);
        return f;
    }

    auto* fBase = new TF1(baseName, func, xmin, xmax, 3);
    auto* fKernel = new TF1(kernelName, GausKernel, -1., 1., 1);

    TF1Convolution conv(fBase, fKernel, xmin, xmax, true);
    conv.SetNofPointsFFT(1000);

    auto* f = new TF1(fullName, conv, xmin, xmax, conv.GetNpar());
    f->SetNpx(kNpxSmooth);
    return f;
}

void RhoFitPanel::RebuildComponents(bool convolved, double xmin, double xmax) {
    delete fRhoFunc;
    delete fF2Func;
    delete fF0Func;
    fRhoFunc = BuildComponentFunction(0, convolved, xmin, xmax);
    fF2Func = BuildComponentFunction(1, convolved, xmin, xmax);
    fF0Func = BuildComponentFunction(2, convolved, xmin, xmax);
}

void RhoFitPanel::SyncComponentParams() {
    if (!fFitFunc) {
        return;
    }

    bool convolved = fConvolveCheck->IsOn();
    double sigma = convolved ? fFitFunc->GetParameter(kNParBase) : 0.;

    auto sync = [&](TF1* f, int baseParIdx) {
        if (!f) {
            return;
        }
        f->SetParameter(0, fFitFunc->GetParameter(baseParIdx));
        f->SetParameter(1, fFitFunc->GetParameter(baseParIdx + 1));
        f->SetParameter(2, fFitFunc->GetParameter(baseParIdx + 2));
        if (convolved) {
            f->SetParameter(3, sigma);
        }
    };

    sync(fRhoFunc, 0);
    sync(fF2Func, 3);
    sync(fF0Func, 6);
}

void RhoFitPanel::ApplyWidgetsToFunction(TF1* f) const {
    bool convolved = fConvolveCheck->IsOn();
    int nPar = convolved ? kNParBase + 1 : kNParBase;

    for (int i = 0; i < nPar; i++) {
        const ParamRow& pr = fParamRows[i];
        double val = pr.value->GetNumber();

        if (pr.fix->IsOn()) {
            f->SetParameter(i, val);
            f->FixParameter(i, val);
        } else {
            double pmin = pr.min->GetNumber();
            double pmax = pr.max->GetNumber();

            if (pmax > pmin) {
                // Keep the starting value strictly inside (pmin, pmax).
                // A bounded Minuit fit's internal parameter transform has
                // an infinite derivative exactly at the boundary, so
                // starting there (which happens routinely here, since the
                // "value" field gets overwritten with the previous fit's
                // result -- itself often sitting right on a limit -- after
                // every fit) is numerically unstable and can let the next
                // fit's result end up outside the declared range entirely.
                double margin = 1.e-6 * (pmax - pmin);
                if (val <= pmin) {
                    val = pmin + margin;
                } else if (val >= pmax) {
                    val = pmax - margin;
                }

                f->SetParameter(i, val);
                f->SetParLimits(i, pmin, pmax);
            } else {
                f->SetParameter(i, val);
            }
        }
    }
}

void RhoFitPanel::Redraw() {
    if (!fHist) {
        return;
    }

    TCanvas* c = fEmbCanvas->GetCanvas();
    c->cd();
    fHist->SetLineColor(kBlack);
    fHist->SetMarkerStyle(20);
    fHist->SetMarkerSize(0.6);
    fHist->Draw("PE");

    delete fLegend;
    fLegend = new TLegend(0.65, 0.65, 0.89, 0.89);
    fLegend->SetBorderSize(0);
    fLegend->SetFillStyle(0);
    fLegend->AddEntry(fHist, fHist->GetName(), "PE");

    if (fFitFunc) {
        fFitFunc->SetNpx(kNpxSmooth);
        fFitFunc->SetLineColor(kBlue);
        fFitFunc->SetLineWidth(2);
        fFitFunc->Draw("same");
        fLegend->AddEntry(fFitFunc, "Fit (#rho + f_{2} + f_{0})", "L");
    }

    if (fRhoFunc) {
        fRhoFunc->SetNpx(kNpxSmooth);
        fRhoFunc->SetLineColor(kRed);
        fRhoFunc->SetLineWidth(2);
        fRhoFunc->Draw("same");
        fLegend->AddEntry(fRhoFunc, "#rho(770)", "L");
    }

    if (fF2Func) {
        fF2Func->SetNpx(kNpxSmooth);
        fF2Func->SetLineColor(kGreen);
        fF2Func->SetLineWidth(2);
        fF2Func->Draw("same");
        fLegend->AddEntry(fF2Func, "f_{2}(1270)", "L");
    }

    if (fF0Func) {
        fF0Func->SetNpx(kNpxSmooth);
        fF0Func->SetLineColor(kMagenta + 1);
        fF0Func->SetLineWidth(2);
        fF0Func->Draw("same");
        fLegend->AddEntry(fF0Func, "f_{0}(1500)", "L");
    }

    fLegend->Draw();

    c->Modified();
    c->Update();
}

void RhoFitPanel::OnHistSelected(Int_t id) {
    if (id < 0 || id >= static_cast<int>(fHistNames.size())) {
        return;
    }
    LoadHistogram(fHistNames[id].c_str());
}

void RhoFitPanel::OnConvolveToggled() {
    bool on = fConvolveCheck->IsOn();
    fParamRows[kNParBase].value->SetState(on);
    fParamRows[kNParBase].min->SetState(on);
    fParamRows[kNParBase].max->SetState(on);
    fParamRows[kNParBase].fix->SetState(on ? kButtonUp : kButtonDisabled);
}

void RhoFitPanel::DoFit() {
    if (!fHist) {
        fStatusLabel->SetText("No histogram selected.");
        return;
    }

    double fitMin = fFitMinEntry->GetNumber();
    double fitMax = fFitMaxEntry->GetNumber();
    double drawMin = fDrawMinEntry->GetNumber();
    double drawMax = fDrawMaxEntry->GetNumber();

    bool convolved = fConvolveCheck->IsOn();

    // Build over the *wide* draw range from the start (matching how the
    // components are built) so that, when convolved, TF1Convolution's FFT
    // grid covers the full domain we will ever evaluate/draw over. FFT
    // convolution is periodic over its construction range; narrowing to the
    // fit range here and widening back afterward via SetRange() would only
    // change what range the outer TF1 reports -- the FFT grid underneath
    // stays fixed at construction time, and evaluating outside it wraps
    // around and picks up the function's behavior from the other edge of
    // the range (visible as the rho threshold's sharp rise "leaking" into
    // the high-mass tail).
    delete fFitFunc;
    fFitFunc = BuildFunction(convolved, drawMin, drawMax);
    ApplyWidgetsToFunction(fFitFunc);

    // Restrict the actual fit to [fitMin, fitMax]: safe to change on the
    // outer TF1 without touching the underlying FFT grid.
    fFitFunc->SetRange(fitMin, fitMax);
    TFitResultPtr result = fHist->Fit(fFitFunc, "RS0");

    int nPar = convolved ? kNParBase + 1 : kNParBase;
    for (int i = 0; i < nPar; i++) {
        fParamRows[i].value->SetNumber(fFitFunc->GetParameter(i));
    }

    // Widen back for display -- safe now, since the FFT grid already covers
    // this range.
    fFitFunc->SetRange(drawMin, drawMax);
    RebuildComponents(convolved, drawMin, drawMax);
    SyncComponentParams();

    // Flag it plainly if any free parameter's fitted value ended up outside
    // the range the user declared for it -- this can still happen for a
    // genuinely difficult/degenerate fit even with a well-behaved starting
    // point, and should never be silently plotted as if nothing were wrong.
    std::string outOfRangeList;
    for (int i = 0; i < nPar; i++) {
        const ParamRow& pr = fParamRows[i];
        if (pr.fix->IsOn()) {
            continue;
        }
        double pmin = pr.min->GetNumber();
        double pmax = pr.max->GetNumber();
        if (pmax <= pmin) {
            continue; // no range was set for this parameter
        }
        double v = fFitFunc->GetParameter(i);
        if (v < pmin || v > pmax) {
            if (!outOfRangeList.empty()) {
                outOfRangeList += ", ";
            }
            outOfRangeList += kParNames[i];
        }
    }

    char buf[512];
    if (result.Get() && result->IsValid()) {
        if (!outOfRangeList.empty()) {
            snprintf(buf, sizeof(buf), "Fit OK (chi2/ndf=%.2f/%d) but OUTSIDE declared range: %s",
                    fFitFunc->GetChisquare(), fFitFunc->GetNDF(), outOfRangeList.c_str());
        } else {
            snprintf(buf, sizeof(buf), "Fit OK: chi2/ndf = %.2f/%d", fFitFunc->GetChisquare(), fFitFunc->GetNDF());
        }
    } else {
        snprintf(buf, sizeof(buf), "Fit did not converge.");
    }
    fStatusLabel->SetText(buf);

    Redraw();
}

void RhoFitPanel::DoDrawFunction() {
    if (!fHist) {
        fStatusLabel->SetText("No histogram selected.");
        return;
    }

    double drawMin = fDrawMinEntry->GetNumber();
    double drawMax = fDrawMaxEntry->GetNumber();
    bool convolved = fConvolveCheck->IsOn();

    delete fFitFunc;
    fFitFunc = BuildFunction(convolved, drawMin, drawMax);
    ApplyWidgetsToFunction(fFitFunc);
    RebuildComponents(convolved, drawMin, drawMax);
    SyncComponentParams();

    fStatusLabel->SetText("Showing function with current parameter values (not fitted).");

    Redraw();
}

void RhoFitPanel::DoSave() {
    if (!fHist) {
        fStatusLabel->SetText("No histogram selected.");
        return;
    }

    gSystem->Exec("mkdir -p Figs");

    std::string base = "Figs/RhoFit_" + std::string(fHist->GetName());

    TCanvas* c = fEmbCanvas->GetCanvas();
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((base + ".root").c_str());

    fStatusLabel->SetText(Form("Saved %s.{pdf,png,root}", base.c_str()));
}

void RhoFitPanel::DoExit() {
    CloseWindow();
    gApplication->Terminate(0);
}

// ================================ Entry point ================================

void RhoFitGUI(const char* fileName = "Hists/Hists_PionPairSkim_All.root") {
    new RhoFitPanel(gClient->GetRoot(), 1250, 720, fileName);
}
