#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>

struct CrystalHistograms {
  TH1F *hist;
  TH1F *zoomedHist;
  TH1F *peakHist;
};

struct Histograms {
  TH1F *hist;
  TH1F *zoomedHist;
  TH1F *peakHist;
  std::vector<CrystalHistograms> crystalHists;
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
    file->Close();
    delete file;
    return {};
  }

  result.hist = static_cast<TH1F *>(hist->Clone());
  result.zoomedHist = static_cast<TH1F *>(zoomedHist->Clone());
  result.peakHist = static_cast<TH1F *>(peakHist->Clone());

  result.hist->SetDirectory(0);
  result.zoomedHist->SetDirectory(0);
  result.peakHist->SetDirectory(0);

  std::cout << filename << ": loaded histograms" << std::endl;

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
  if (!Constants::NORMALIZE_BY_TIME) {
    std::cerr << "ERROR: SubtractBackground requires NORMALIZE_BY_TIME = kTRUE"
              << std::endl;
    return;
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  std::vector<DatasetPair> datasets;

  datasets.push_back({"CdShield_10_01132026",
                      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CdShield_25_01132026",
                      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113});

  datasets.push_back({"CuShield_10_01132026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CuShield_10_01142026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114});

  datasets.push_back({"NoShield_5_01152026",
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115});

  datasets.push_back(
      {"NoShield_GraphiteCastle_10_01162026",
       Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
       Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116});

  const Int_t N_CRYSTALS = 4;
  std::vector<Histograms> bkgSubtractedResults;

  for (Int_t i = 0; i < (Int_t)datasets.size(); i++) {
    std::cout << "Processing: " << datasets[i].name << std::endl;

    Histograms signal = GetHistograms(datasets[i].signalFile);
    Histograms background = GetHistograms(datasets[i].backgroundFile);

    Histograms bkgSubtracted =
        SubtractHistograms(signal, background, datasets[i].name);

    TFile *sigFile =
        new TFile("root_files/" + datasets[i].signalFile + ".root", "READ");
    TFile *bkgFile =
        new TFile("root_files/" + datasets[i].backgroundFile + ".root", "READ");

    for (Int_t c = 0; c < N_CRYSTALS; c++) {
      TString histName = Form("hist_crystal%d", c);
      TString zoomedName = Form("zoomedHist_crystal%d", c);
      TString peakName = Form("peakHist_crystal%d", c);

      TH1F *sigHist = static_cast<TH1F *>(sigFile->Get(histName));
      TH1F *sigZoomed = static_cast<TH1F *>(sigFile->Get(zoomedName));
      TH1F *sigPeak = static_cast<TH1F *>(sigFile->Get(peakName));
      TH1F *bkgHist = static_cast<TH1F *>(bkgFile->Get(histName));
      TH1F *bkgZoomed = static_cast<TH1F *>(bkgFile->Get(zoomedName));
      TH1F *bkgPeak = static_cast<TH1F *>(bkgFile->Get(peakName));

      if (!sigHist || !bkgHist) {
        std::cerr << "WARNING: missing crystal " << c << " histograms for "
                  << datasets[i].name << std::endl;
        continue;
      }

      TH1F *subHist = static_cast<TH1F *>(sigHist->Clone(
          Form("hist_%s_crystal%d_BkgSubtracted", datasets[i].name.Data(), c)));
      subHist->Add(bkgHist, -1);
      subHist->SetDirectory(0);

      TH1F *subZoomed = static_cast<TH1F *>(
          sigZoomed->Clone(Form("zoomedHist_%s_crystal%d_BkgSubtracted",
                                datasets[i].name.Data(), c)));
      subZoomed->Add(bkgZoomed, -1);
      subZoomed->SetDirectory(0);

      TH1F *subPeak = static_cast<TH1F *>(sigPeak->Clone(Form(
          "peakHist_%s_crystal%d_BkgSubtracted", datasets[i].name.Data(), c)));
      subPeak->Add(bkgPeak, -1);
      subPeak->SetDirectory(0);

      CrystalHistograms ch = {subHist, subZoomed, subPeak};
      bkgSubtracted.crystalHists.push_back(ch);
    }

    sigFile->Close();
    bkgFile->Close();
    delete sigFile;
    delete bkgFile;

    bkgSubtractedResults.push_back(bkgSubtracted);
  }

  TH1F *hist_Combined =
      static_cast<TH1F *>(bkgSubtractedResults[0].hist->Clone("hist_Combined"));
  hist_Combined->Reset();

  TH1F *zoomedHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].zoomedHist->Clone("zoomedHist_Combined"));
  zoomedHist_Combined->Reset();

  TH1F *peakHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].peakHist->Clone("peakHist_Combined"));
  peakHist_Combined->Reset();

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    hist_Combined->Add(bkgSubtractedResults[i].hist);
    zoomedHist_Combined->Add(bkgSubtractedResults[i].zoomedHist);
    peakHist_Combined->Add(bkgSubtractedResults[i].peakHist);
  }

  hist_Combined->SetTitle(
      Form("Background Subtracted; Energy [keV]; Counts / %d eV / s",
           Constants::BIN_WIDTH_EV));

  zoomedHist_Combined->SetTitle(
      Form("Background Subtracted (Zoomed); Energy [keV]; Counts / %d eV / s",
           Constants::BIN_WIDTH_EV));

  peakHist_Combined->SetTitle(Form(
      "Background Subtracted (Peak Region); Energy [keV]; Counts / %d eV / s",
      Constants::BIN_WIDTH_EV));

  TCanvas *canvasFull = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->SetFillStyle(0);
  hist_Combined->SetLineWidth(2);
  hist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasFull, "background_subtracted",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TCanvas *canvasZoomed = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->SetFillStyle(0);
  zoomedHist_Combined->SetLineWidth(2);
  zoomedHist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasZoomed, "background_subtracted_zoomed",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TCanvas *canvasPeak = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(peakHist_Combined, kP10Violet);
  peakHist_Combined->SetFillStyle(0);
  peakHist_Combined->SetLineWidth(2);
  peakHist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasPeak, "background_subtracted_peak",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TString outputName = "Combined_BkgSubtracted_Normalized";
  TFile *outputFile =
      new TFile("root_files/" + outputName + ".root", "RECREATE");

  hist_Combined->Write("histCombined", TObject::kOverwrite);
  zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
  peakHist_Combined->Write("peakHistCombined", TObject::kOverwrite);

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    bkgSubtractedResults[i].hist->Write();
    bkgSubtractedResults[i].zoomedHist->Write();
    bkgSubtractedResults[i].peakHist->Write();

    for (Int_t c = 0; c < (Int_t)bkgSubtractedResults[i].crystalHists.size();
         c++) {
      bkgSubtractedResults[i].crystalHists[c].hist->Write();
      bkgSubtractedResults[i].crystalHists[c].zoomedHist->Write();
      bkgSubtractedResults[i].crystalHists[c].peakHist->Write();
    }
  }

  outputFile->Close();
  delete outputFile;

  std::cout << "Background subtraction and combination complete!" << std::endl;
  std::cout << "Output saved to: root_files/" << outputName << ".root"
            << std::endl;
}
