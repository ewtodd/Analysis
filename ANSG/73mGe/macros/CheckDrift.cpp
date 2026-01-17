#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <vector>

void CheckDrift() {
  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  std::vector<TString> calibFiles = {
      "01132026-PostReactor-Calibration", "01152026-NewSetup-PostReactor-Am241",
      "01152026-NewSetup-PostReactor-Ba133",
      "01162026-NoShield-PostReactor-Am241-Ba133"};

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  TCanvas *c1 = new TCanvas("c1", "Calibration Spectra Comparison", 1200, 800);

  TLegend *legend = new TLegend(0.42, 0.60, 0.89, 0.89);
  legend->SetBorderSize(1);
  legend->SetTextSize(0.03);

  std::vector<TH1F *> histograms;
  Double_t maxVal = 0;

  for (size_t i = 0; i < calibFiles.size(); i++) {
    TString filepath = "root_files/" + calibFiles[i] + suffix + ".root";
    TFile *file = TFile::Open(filepath, "READ");

    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: Could not open file " << filepath << std::endl;
      continue;
    }

    TH1F *hist = (TH1F *)file->Get("zoomedHist");

    if (!hist) {
      std::cerr << "ERROR: Could not find zoomedHist in " << filepath
                << std::endl;
      file->Close();
      continue;
    }

    TH1F *histClone = (TH1F *)hist->Clone(Form("hist_%zu", i));
    histClone->SetDirectory(0);

    Double_t integral = histClone->Integral();
    if (integral > 0) {
      histClone->Scale(1.0 / integral);
    }

    histClone->SetLineColor(colors[i % colors.size()]);
    histClone->SetLineWidth(2);
    histClone->SetFillStyle(0);
    histClone->SetStats(0);

    histClone->GetYaxis()->SetTitle(
        Form("Normalized Counts / %d eV", Constants::BIN_WIDTH_EV));
    histClone->GetYaxis()->SetTitleOffset(1.3);

    histograms.push_back(histClone);

    Double_t histMax = histClone->GetMaximum();
    if (histMax > maxVal)
      maxVal = histMax;

    legend->AddEntry(histClone, calibFiles[i], "l");

    file->Close();
  }

  if (histograms.size() > 0) {
    histograms[0]->SetMaximum(maxVal * 1.5);
    histograms[0]->Draw("HIST");

    for (size_t i = 1; i < histograms.size(); i++) {
      histograms[i]->Draw("HIST SAME");
    }

    legend->Draw();

    c1->Update();
    PlottingUtils::SaveFigure(c1, "drift_check.png", kFALSE);

    std::cout << "Plotted " << histograms.size() << " calibration spectra"
              << std::endl;
  } else {
    std::cerr << "ERROR: No histograms were loaded successfully" << std::endl;
  }
}
