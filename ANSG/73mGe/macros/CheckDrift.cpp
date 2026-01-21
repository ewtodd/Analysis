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

  std::vector<TString> calibFiles = {
      Constants::POSTREACTOR_AM241_01132026,
      Constants::POSTREACTOR_AM241_01152026,
      Constants::POSTREACTOR_BA133_01152026,
      Constants::POSTREACTOR_AM241_BA133_01162026};

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  TCanvas *c1 = new TCanvas("c1", "Calibration Spectra Comparison", 1200, 800);
  PlottingUtils::ConfigureCanvas(c1, kFALSE);

  TLegend *legend = new TLegend(0.42, 0.60, 0.89, 0.89);
  legend->SetBorderSize(1);
  legend->SetTextSize(0.03);

  std::vector<TH1F *> histograms;
  Double_t maxVal = 0;

  for (size_t i = 0; i < calibFiles.size(); i++) {
    TString filepath = "root_files/" + calibFiles[i] + ".root";
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
    histograms[0]->SetMaximum(maxVal * 1.2);
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
