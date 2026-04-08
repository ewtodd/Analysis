#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <cmath>
#include <iomanip>
#include <vector>

TH1D *LoadBkgSubtractedHistogram(const TString filename,
                                 const TString histName) {
  TFile *file = new TFile("root_files/" + filename + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filename << ".root" << std::endl;
    return nullptr;
  }

  TH1D *hist = static_cast<TH1D *>(file->Get(histName));
  if (!hist) {
    std::cerr << "ERROR: Cannot find " << histName << " in " << filename
              << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }
  hist->SetDirectory(0);
  file->Close();
  delete file;
  return hist;
}

FitResult FitGePeak(TH1D *hist, const TString label, const Bool_t interactive) {
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter =
      new FittingUtils(hist, 64, 77, use_flat_background, use_step,
                       use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter->SetInteractive();

  FitResult result = fitter->FitSinglePeak(label, "Ge_68.75keV");
  delete fitter;
  return result;
}

void FitsBkgSubtracted() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  TString outputFile = "Combined_BkgSubtracted_Normalized";

  std::vector<TString> datasetNames = {
      //  "CdShield_10_01132026",  "CdShield_25_01132026",
      //  "CuShield_10_01132026",  "CuShield_10_01142026",
      "NoShield_5_01152026", "NoShield_GraphiteCastle_10_01162026"};

  std::vector<TString> run_names;
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> reduced_chi2;

  for (Int_t i = 0; i < (Int_t)datasetNames.size(); i++) {
    TString histName =
        Form("zoomedHist_%s_BkgSubtracted", datasetNames[i].Data());
    TH1D *hist = LoadBkgSubtractedHistogram(outputFile, histName);
    if (!hist)
      continue;

    FitResult result = FitGePeak(hist, datasetNames[i], interactive);
    delete hist;

    if (result.peaks.empty())
      continue;

    run_names.push_back(datasetNames[i]);
    mu.push_back(result.peaks.at(0).mu);
    mu_errors.push_back(result.peaks.at(0).mu_error);
    reduced_chi2.push_back(result.reduced_chi2);
  }

  std::cout << std::endl;
  std::cout << "Background-Subtracted Ge Peak Results:" << std::endl;

  for (Int_t i = 0; i < (Int_t)mu.size(); i++) {
    std::cout << std::left << std::setw(50) << run_names[i] << ": "
              << std::fixed << std::setprecision(4) << mu[i] << " +/- "
              << mu_errors[i] << " keV"
              << " (chi2/ndf = " << std::setprecision(3) << reduced_chi2[i]
              << ")" << std::endl;
  }

  Float_t sum_weights = 0.0;
  Float_t weighted_sum = 0.0;
  for (Int_t i = 0; i < (Int_t)mu.size(); i++) {
    if (mu_errors[i] > 0) {
      Float_t weight = 1.0 / (mu_errors[i] * mu_errors[i]);
      weighted_sum += mu[i] * weight;
      sum_weights += weight;
    }
  }
  if (sum_weights > 0) {
    Float_t combined_mu = weighted_sum / sum_weights;
    Float_t combined_error = std::sqrt(1.0 / sum_weights);
    std::cout << std::endl;
    std::cout << "Combined Ge mu (bkg-subtracted): " << std::fixed
              << std::setprecision(4) << combined_mu << " +/- "
              << combined_error << " keV" << std::endl;
  }
}
