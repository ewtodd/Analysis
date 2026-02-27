#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

struct Histograms {
  TH1F *hist;
  TH1F *zoomedHist;
  TH1F *peakHist;
};

struct DatasetPair {
  TString name;
  TString signalFile;
  TString backgroundFile;
};

Histograms GetHistograms(TString filename) {
  Histograms result;

  TFile *file = new TFile("root_files/" + filename + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filename << ".root" << std::endl;
    return {};
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("hist"));
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get("zoomedHist"));
  TH1F *peakHist = static_cast<TH1F *>(file->Get("peakHist"));

  if (!hist || !zoomedHist || !peakHist) {
    std::cerr << "ERROR: Cannot retrieve histograms from " << filename
              << std::endl;
    return {};
  }

  result.hist = static_cast<TH1F *>(hist->Clone());
  result.zoomedHist = static_cast<TH1F *>(zoomedHist->Clone());
  result.peakHist = static_cast<TH1F *>(peakHist->Clone());

  result.hist->SetDirectory(0);
  result.zoomedHist->SetDirectory(0);
  result.peakHist->SetDirectory(0);

  file->Close();
  delete file;

  return result;
}

Histograms SubtractHistograms(const Histograms &signal,
                              const Histograms &background, TString name) {
  Histograms result;

  result.hist = static_cast<TH1F *>(
      signal.hist->Clone(Form("hist_%s_BkgSubtracted", name.Data())));
  result.hist->Add(background.hist, -1);

  result.zoomedHist = static_cast<TH1F *>(signal.zoomedHist->Clone(
      Form("zoomedHist_%s_BkgSubtracted", name.Data())));
  result.zoomedHist->Add(background.zoomedHist, -1);

  result.peakHist = static_cast<TH1F *>(
      signal.peakHist->Clone(Form("peakHist_%s_BkgSubtracted", name.Data())));
  result.peakHist->Add(background.peakHist, -1);

  return result;
}

void SubtractBackground() {
  InitUtils::SetROOTPreferences();

  if (Constants::NORMALIZE_BY_TIME) {
    std::cout << "Running in NORMALIZED mode: performing background subtraction"
              << std::endl;
    if (Constants::FILTERED) {
      std::cout << "WARNING: normalization by time on filtered data may not be "
                   "reliable."
                << std::endl;
    }
  }

  std::vector<DatasetPair> datasets;

  datasets.push_back({"CdShield_10_20260113",
                      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CdShield_25_20260113",
                      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113});

  datasets.push_back({"CuShield_10_20260113",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CuShield_10_20260114",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114});

  datasets.push_back({"NoShield_5_20260115",
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115});

  datasets.push_back(
      {"NoShield_GraphiteCastle_10_20260116",
       Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
       Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116});

  std::vector<Histograms> bkgSubtractedResults;

  for (const auto &dataset : datasets) {
    std::cout << "Processing: " << dataset.name << std::endl;

    Histograms signal = GetHistograms(dataset.signalFile);
    Histograms background = GetHistograms(dataset.backgroundFile);

    Histograms bkgSubtracted =
        SubtractHistograms(signal, background, dataset.name);
    bkgSubtractedResults.push_back(bkgSubtracted);
  }

  TString perSecond = Constants::NORMALIZE_BY_TIME ? " / s" : "";

  TH1F *hist_Combined =
      static_cast<TH1F *>(bkgSubtractedResults[0].hist->Clone("hist_Combined"));
  hist_Combined->Reset();

  TH1F *zoomedHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].zoomedHist->Clone("zoomedHist_Combined"));
  zoomedHist_Combined->Reset();

  TH1F *peakHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].peakHist->Clone("peakHist_Combined"));
  peakHist_Combined->Reset();

  for (const auto &result : bkgSubtractedResults) {
    hist_Combined->Add(result.hist);
    zoomedHist_Combined->Add(result.zoomedHist);
    peakHist_Combined->Add(result.peakHist);
  }

  hist_Combined->SetTitle(
      Form("Background Subtracted; Energy [keV]; Counts / %d eV%s",
           Constants::BIN_WIDTH_EV, perSecond.Data()));

  zoomedHist_Combined->SetTitle(
      Form("Background Subtracted (Zoomed); Energy [keV]; Counts / %d eV%s",
           Constants::BIN_WIDTH_EV, perSecond.Data()));

  peakHist_Combined->SetTitle(Form(
      "Background Subtracted (Peak Region); Energy [keV]; Counts / %d eV%s",
      Constants::BIN_WIDTH_EV, perSecond.Data()));

  TCanvas *canvasFull = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->Draw("HIST");
  hist_Combined->SetFillStyle(0);
  hist_Combined->SetLineWidth(2);
  PlottingUtils::SaveFigure(canvasFull, "background_subtracted",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvasZoomed = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->Draw("HIST");
  zoomedHist_Combined->SetFillStyle(0);
  zoomedHist_Combined->SetLineWidth(2);
  PlottingUtils::SaveFigure(canvasZoomed, "background_subtracted_zoomed",
                            PlotSaveOptions::kLINEAR);

  TCanvas *canvasPeak = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(peakHist_Combined, kP10Violet);
  peakHist_Combined->Draw("HIST");
  peakHist_Combined->SetFillStyle(0);
  peakHist_Combined->SetLineWidth(2);
  PlottingUtils::SaveFigure(canvasPeak, "background_subtracted_peak",
                            PlotSaveOptions::kLINEAR);

  TString suffix = Constants::FILTERED ? "_Filtered" : "";
  if (Constants::NORMALIZE_BY_TIME) {
    suffix += "_Normalized";
  }

  TString outputName = "Combined_BkgSubtracted" + suffix;
  TFile *outputFile =
      new TFile("root_files/" + outputName + ".root", "RECREATE");

  hist_Combined->Write("histCombined", TObject::kOverwrite);
  zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
  peakHist_Combined->Write("peakHistCombined", TObject::kOverwrite);

  for (size_t i = 0; i < bkgSubtractedResults.size(); i++) {
    bkgSubtractedResults[i].hist->Write();
    bkgSubtractedResults[i].zoomedHist->Write();
    bkgSubtractedResults[i].peakHist->Write();
  }

  outputFile->Close();
  delete outputFile;

  std::cout << "Background subtraction and combination complete!" << std::endl;
  std::cout << "Output saved to: root_files/" << outputName << ".root"
            << std::endl;
}
