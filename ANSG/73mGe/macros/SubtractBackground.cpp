#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "HyperEMGFitHelpers.hpp"
#include "InitUtils.hpp"
#include "InteractiveSubtractionEditor.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <iomanip>

const Float_t RESIDUAL_REGION_LOW = 72.0;
const Float_t RESIDUAL_REGION_HIGH = 90.0;

struct DatasetPair {
  TString name;
  TString signalFile;
  TString backgroundFile;
};

Float_t ComputeTotalLiveTime(TFile *file) {
  TString prefix = Constants::USE_FILTERED ? "Filtered" : "Unfiltered";
  Float_t total = 0;
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TString paramName = Form("LiveTime_%s_Crystal%d_s", prefix.Data(), c);
    TParameter<Float_t> *param =
        static_cast<TParameter<Float_t> *>(file->Get(paramName));
    if (param)
      total += param->GetVal();
  }
  return total;
}

void SubtractBackground() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

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

  struct BkgSubResult {
    TH1F *hist;
    TH1F *zoomedHist;
    TH1F *peakHist;
    TString name;
    HyperEMGFitResult fit;
  };

  std::vector<BkgSubResult> bkgSubtractedResults;

  for (Int_t i = 0; i < (Int_t)datasets.size(); i++) {
    std::cout << "Processing: " << datasets[i].name << std::endl;

    TFile *sigFile =
        new TFile("root_files/" + datasets[i].signalFile + ".root", "READ");
    TFile *bkgFile =
        new TFile("root_files/" + datasets[i].backgroundFile + ".root", "READ");

    if (!sigFile || sigFile->IsZombie() || !bkgFile || bkgFile->IsZombie()) {
      std::cerr << "ERROR: Cannot open files for " << datasets[i].name
                << std::endl;
      if (sigFile)
        delete sigFile;
      if (bkgFile)
        delete bkgFile;
      continue;
    }

    // Total live times for normalization
    Float_t sigLT = ComputeTotalLiveTime(sigFile);
    Float_t bkgLT = ComputeTotalLiveTime(bkgFile);
    Float_t ltRatio = (bkgLT > 0) ? sigLT / bkgLT : 1.0;

    std::cout << "  LT ratio = " << std::fixed << std::setprecision(4)
              << ltRatio << " (sig = " << sigLT << " s, bkg = " << bkgLT
              << " s)" << std::endl;

    // Load calibrated (SNIP-subtracted + calibrated) histograms
    TH1F *sigHist =
        static_cast<TH1F *>(sigFile->Get("calibrated_hist"));
    TH1F *sigZoomed =
        static_cast<TH1F *>(sigFile->Get("calibrated_zoomedHist"));
    TH1F *sigPeak =
        static_cast<TH1F *>(sigFile->Get("calibrated_peakHist"));
    TH1F *bkgHist =
        static_cast<TH1F *>(bkgFile->Get("calibrated_hist"));
    TH1F *bkgZoomed =
        static_cast<TH1F *>(bkgFile->Get("calibrated_zoomedHist"));
    TH1F *bkgPeak =
        static_cast<TH1F *>(bkgFile->Get("calibrated_peakHist"));

    if (!sigHist || !sigZoomed || !sigPeak || !bkgHist || !bkgZoomed ||
        !bkgPeak) {
      std::cerr << "ERROR: Missing calibrated histograms for "
                << datasets[i].name << std::endl;
      sigFile->Close();
      bkgFile->Close();
      delete sigFile;
      delete bkgFile;
      continue;
    }

    sigHist->SetDirectory(0);
    sigZoomed->SetDirectory(0);
    sigPeak->SetDirectory(0);
    bkgHist->SetDirectory(0);
    bkgZoomed->SetDirectory(0);
    bkgPeak->SetDirectory(0);

    // Find optimal correction scale in residual region
    // Normalize background by live time ratio first
    TH1F *bkgZoomedNormed = static_cast<TH1F *>(
        bkgZoomed->Clone(PlottingUtils::GetRandomName()));
    bkgZoomedNormed->SetDirectory(0);
    bkgZoomedNormed->Scale(ltRatio);

    Double_t sb = 0, bb = 0;
    for (Int_t bin = 1; bin <= sigZoomed->GetNbinsX(); bin++) {
      Double_t x = sigZoomed->GetBinCenter(bin);
      if (x >= RESIDUAL_REGION_LOW && x <= RESIDUAL_REGION_HIGH) {
        Double_t s = sigZoomed->GetBinContent(bin);
        Double_t b = bkgZoomedNormed->GetBinContent(bin);
        sb += s * b;
        bb += b * b;
      }
    }
    Double_t best_correction = (bb > 0) ? sb / bb : 1.0;

    std::cout << "  Optimal correction = " << std::fixed << std::setprecision(4)
              << best_correction << std::endl;

    // Interactive editor
    Double_t correction = best_correction;
    if (interactive) {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      Double_t editor_result = LaunchSubtractionEditor(
          sigZoomed, bkgZoomedNormed, best_correction,
          RESIDUAL_REGION_LOW, RESIDUAL_REGION_HIGH, datasets[i].name);
      gROOT->SetBatch(was_batch);
      if (editor_result > 0)
        correction = editor_result;
      std::cout << "  Correction factor: " << std::fixed << std::setprecision(4)
                << correction << std::endl;
    }

    delete bkgZoomedNormed;

    // Final subtraction: signal - correction * ltRatio * background
    Float_t scale = correction * ltRatio;

    BkgSubResult result;

    result.hist = static_cast<TH1F *>(sigHist->Clone(
        Form("hist_%s_BkgSubtracted", datasets[i].name.Data())));
    result.hist->SetDirectory(0);
    result.hist->Add(bkgHist, -scale);
    result.hist->Scale(1.0 / sigLT);

    result.zoomedHist = static_cast<TH1F *>(sigZoomed->Clone(
        Form("zoomedHist_%s_BkgSubtracted", datasets[i].name.Data())));
    result.zoomedHist->SetDirectory(0);
    result.zoomedHist->Add(bkgZoomed, -scale);
    result.zoomedHist->Scale(1.0 / sigLT);

    result.peakHist = static_cast<TH1F *>(sigPeak->Clone(
        Form("peakHist_%s_BkgSubtracted", datasets[i].name.Data())));
    result.peakHist->SetDirectory(0);
    result.peakHist->Add(bkgPeak, -scale);
    result.peakHist->Scale(1.0 / sigLT);

    sigFile->Close();
    bkgFile->Close();
    delete sigFile;
    delete bkgFile;

    TCanvas *canvasPair = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureHistogram(result.peakHist, kP10Violet);
    result.peakHist->Sumw2(0);
    result.peakHist->SetMarkerStyle(20);
    result.peakHist->SetMarkerSize(0.65);
    result.peakHist->Draw("P");
    PlottingUtils::SaveFigure(canvasPair, "peak_" + datasets[i].name,
                              "backgroundSubtraction/pairs",
                              PlotSaveOptions::kLINEAR);

    // Fit Ge-73m peak in this dataset
    result.name = datasets[i].name;
    result.fit = FitSingleHyperEMG(result.zoomedHist, 64.0, 77.0, kTRUE,
                                    datasets[i].name, "Ge73m_68keV",
                                    interactive);
    if (result.fit.valid) {
      std::cout << "  " << datasets[i].name << ": mu = " << std::fixed
                << std::setprecision(4) << result.fit.peaks[0].mu << " +/- "
                << result.fit.peaks[0].mu_error
                << " keV, sigma = " << result.fit.peaks[0].sigma << " +/- "
                << result.fit.peaks[0].sigma_error
                << " keV, chi2/ndf = " << std::setprecision(3)
                << result.fit.reduced_chi2 << std::endl;
    } else {
      std::cout << "  " << datasets[i].name << ": FIT FAILED" << std::endl;
    }

    bkgSubtractedResults.push_back(result);
  }

  // Combine all datasets
  TH1F *hist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].hist->Clone("hist_Combined"));
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
      Form("Background Subtracted; Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));
  zoomedHist_Combined->SetTitle(
      Form("Background Subtracted (Zoomed); Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));
  peakHist_Combined->SetTitle(
      Form("Background Subtracted (Peak Region); Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));

  TCanvas *canvasFull = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->SetMarkerStyle(20);
  hist_Combined->SetMarkerSize(0.65);
  hist_Combined->SetMarkerColor(kP10Violet);
  hist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasFull, "background_subtracted",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  TCanvas *canvasZoomed = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->SetMarkerStyle(20);
  zoomedHist_Combined->SetMarkerSize(0.65);
  zoomedHist_Combined->SetMarkerColor(kP10Violet);
  zoomedHist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasZoomed, "background_subtracted_zoomed",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  TCanvas *canvasPeak = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(peakHist_Combined, kP10Violet);
  peakHist_Combined->Sumw2(0);
  peakHist_Combined->SetMarkerStyle(20);
  peakHist_Combined->SetMarkerSize(0.65);
  peakHist_Combined->SetMarkerColor(kP10Violet);
  peakHist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasPeak, "background_subtracted_peak",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  // Fit Ge-73m 68.75 keV peak in the combined spectrum
  std::cout << std::endl;
  std::cout
      << "  Fitting Ge-73m peak in combined background-subtracted spectrum:"
      << std::endl;

  HyperEMGFitResult ge_fit =
      FitSingleHyperEMG(zoomedHist_Combined, 64.0, 77.0, kTRUE, "Combined",
                          "Ge73m_68keV", interactive);
  if (ge_fit.valid) {
    std::cout << "    mu = " << std::fixed << std::setprecision(4)
              << ge_fit.peaks[0].mu << " +/- " << ge_fit.peaks[0].mu_error
              << " keV, sigma = " << ge_fit.peaks[0].sigma << " +/- "
              << ge_fit.peaks[0].sigma_error
              << " keV, chi2/ndf = " << std::setprecision(3)
              << ge_fit.reduced_chi2 << std::endl;
  } else {
    std::cout << "    FIT FAILED" << std::endl;
  }

  // Uncertainty-weighted average from per-dataset fits
  Double_t sum_w = 0, sum_wmu = 0;
  Double_t sum_wsigma = 0;
  Int_t n_valid = 0;

  std::cout << std::endl;
  std::cout << "  Per-dataset Ge-73m fit summary:" << std::endl;
  std::cout << "  " << std::left << std::setw(40) << "Dataset" << std::right
            << std::setw(12) << "mu [keV]" << std::setw(12) << "error"
            << std::setw(12) << "sigma" << std::setw(10) << "chi2/ndf"
            << std::endl;

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    FitResult &fit = bkgSubtractedResults[i].fit;
    if (fit.valid) {
      Double_t mu = fit.peaks[0].mu;
      Double_t mu_err = fit.peaks[0].mu_error;
      Double_t sigma = fit.peaks[0].sigma;
      Double_t w = 1.0 / (mu_err * mu_err);

      sum_w += w;
      sum_wmu += w * mu;
      sum_wsigma += w * sigma;
      n_valid++;

      std::cout << "  " << std::left << std::setw(40)
                << bkgSubtractedResults[i].name << std::right << std::fixed
                << std::setprecision(4) << std::setw(12) << mu << std::setw(12)
                << mu_err << std::setw(12) << sigma << std::setprecision(3)
                << std::setw(10) << fit.reduced_chi2 << std::endl;
    } else {
      std::cout << "  " << std::left << std::setw(40)
                << bkgSubtractedResults[i].name << "  FAILED" << std::endl;
    }
  }

  if (n_valid > 0) {
    Double_t weighted_mu = sum_wmu / sum_w;
    Double_t weighted_mu_err = TMath::Sqrt(1.0 / sum_w);
    Double_t weighted_sigma = sum_wsigma / sum_w;

    std::cout << std::endl;
    std::cout << "  Weighted average: mu = " << std::fixed
              << std::setprecision(4) << weighted_mu << " +/- "
              << weighted_mu_err << " keV, sigma = " << weighted_sigma
              << " keV (" << n_valid << " datasets)" << std::endl;

    if (ge_fit.valid) {
      Double_t combined_mu = ge_fit.peaks[0].mu;
      Double_t combined_err = ge_fit.peaks[0].mu_error;
      Double_t diff = combined_mu - weighted_mu;
      Double_t diff_sigma =
          diff / TMath::Sqrt(combined_err * combined_err +
                             weighted_mu_err * weighted_mu_err);

      std::cout << "  Combined fit:     mu = " << combined_mu << " +/- "
                << combined_err << " keV" << std::endl;
      std::cout << "  Difference:       " << std::setprecision(4) << diff
                << " keV (" << std::setprecision(2) << diff_sigma << " sigma)"
                << std::endl;
    }
  }
  std::cout << std::endl;

  TString outputName = "Combined_BkgSubtracted";
  TFile *outputFile =
      new TFile("root_files/" + outputName + ".root", "RECREATE");

  hist_Combined->Write("histCombined", TObject::kOverwrite);
  zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
  peakHist_Combined->Write("peakHistCombined", TObject::kOverwrite);

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
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
