#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iomanip>
#include <vector>

FitResult FitSinglePeak(const TString input_name, const TString peak_name,
                        const Float_t expected_mu,
                        const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResult result;

  if (peak_name == "Background") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115) {
      fitter = new FittingUtils(zoomedHist, 67, 77, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name ==
               Constants::
                   NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116) {
      fitter = new FittingUtils(zoomedHist, 67, 80, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    }
  }

  result = fitter->FitPeak(input_name, peak_name);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                        const Float_t mu1_init, const Float_t mu2_init,
                        const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResult result;

  if (peak_name == "Pb_KAlpha") {
    if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114) {
      fitter = new FittingUtils(zoomedHist, 60, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CDSHIELDBACKGROUND_10PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 63, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result = fitter->FitDoublePeak(input_name, peak_name, mu1_init, mu2_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResult FitDoublePeakConstrained(const TString input_name,
                                   const TString peak_name,
                                   const PeakFitResult &constrained_peak,
                                   const Float_t mu2_init,
                                   const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResult result;

  if (peak_name == "Ge") {
    if (input_name == Constants::NOSHIELDSIGNAL_5PERCENT_20260115) {
      fitter = new FittingUtils(zoomedHist, 64, 77, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name ==
               Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116) {
      fitter = new FittingUtils(zoomedHist, 60, 77, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    }
  }

  result =
      fitter->FitDoublePeak(input_name, peak_name, constrained_peak, mu2_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                        const FitResult &constrained_peaks,
                        const Float_t mu3_init,
                        const Bool_t use_calibrated = kFALSE) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *zoomedHist = static_cast<TH1F *>(file->Get(histName));

  zoomedHist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResult result;

  if (peak_name == "Ge") {
    if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CDSHIELDSIGNAL_10PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 64, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_20260114) {
      fitter = new FittingUtils(zoomedHist, 62, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else if (input_name == Constants::CUSHIELDSIGNAL_90PERCENT_20260114) {
      fitter = new FittingUtils(zoomedHist, 63, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result =
      fitter->FitTriplePeak(input_name, peak_name, constrained_peaks, mu3_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

void Results() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;

  FitResult calibrated_cd_shield_background_10_percent_20260113 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694, kTRUE);

  FitResult calibrated_cd_shield_background_25_percent_20260113 =
      FitDoublePeak(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694, kTRUE);

  FitResult calibrated_cu_shield_background_10_percent_20260113 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                    "Pb_KAlpha", 72.8042, 74.9694, kTRUE);

  FitResult calibrated_cu_shield_background_10_percent_20260114 =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                    "Pb_KAlpha", 72.8042, 74.9694, kTRUE);

  FitResult calibrated_cd_shield_signal_10_percent_20260113 = FitTriplePeak(
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113, "Ge",
      calibrated_cd_shield_background_10_percent_20260113, 68.75, kTRUE);

  FitResult calibrated_cd_shield_signal_25_percent_20260113 = FitTriplePeak(
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113, "Ge",
      calibrated_cd_shield_background_25_percent_20260113, 68.75, kTRUE);

  FitResult calibrated_cu_shield_signal_10_percent_20260113 = FitTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113, "Ge",
      calibrated_cu_shield_background_10_percent_20260113, 68.75, kTRUE);

  FitResult calibrated_cu_shield_signal_10_percent_20260114 = FitTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, "Ge",
      calibrated_cu_shield_background_10_percent_20260114, 68.75, kTRUE);

  FitResult calibrated_cu_shield_signal_90_percent_20260114 = FitTriplePeak(
      Constants::CUSHIELDSIGNAL_90PERCENT_20260114, "Ge",
      calibrated_cu_shield_background_10_percent_20260114, 68.75, kTRUE);

  FitResult calibrated_no_shield_background_5_percent_20260115 = FitSinglePeak(
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115, "Background", 72, kTRUE);

  FitResult
      calibrated_no_shield_graphite_castle_background_10_percent_20260116 =
          FitSinglePeak(
              Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
              "Background", 72, kTRUE);

  FitResult calibrated_no_shield_signal_5_percent_20260115 =
      FitDoublePeakConstrained(
          Constants::NOSHIELDSIGNAL_5PERCENT_20260115, "Ge",
          calibrated_no_shield_background_5_percent_20260115.peaks.at(0), 68.75,
          kTRUE);

  FitResult calibrated_no_shield_graphite_castle_signal_10_percent_20260116 =
      FitDoublePeakConstrained(
          Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116, "Ge",
          calibrated_no_shield_graphite_castle_background_10_percent_20260116
              .peaks.at(0),
          68.75, kTRUE);

  run_names.push_back("Cd Shield Signal 10% (01/13/2026)");
  mu.push_back(calibrated_cd_shield_signal_10_percent_20260113.peaks.at(2).mu);
  mu_errors.push_back(
      calibrated_cd_shield_signal_10_percent_20260113.peaks.at(2).mu_error);
  reduced_chi2.push_back(
      calibrated_cd_shield_signal_10_percent_20260113.reduced_chi2);

  run_names.push_back("Cd Shield Signal 25% (01/13/2026)");
  mu.push_back(calibrated_cd_shield_signal_25_percent_20260113.peaks.at(2).mu);
  mu_errors.push_back(
      calibrated_cd_shield_signal_25_percent_20260113.peaks.at(2).mu_error);
  reduced_chi2.push_back(
      calibrated_cd_shield_signal_25_percent_20260113.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% (01/13/2026)");
  mu.push_back(calibrated_cu_shield_signal_10_percent_20260113.peaks.at(2).mu);
  mu_errors.push_back(
      calibrated_cu_shield_signal_10_percent_20260113.peaks.at(2).mu_error);
  reduced_chi2.push_back(
      calibrated_cu_shield_signal_10_percent_20260113.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% (01/14/2026)");
  mu.push_back(calibrated_cu_shield_signal_10_percent_20260114.peaks.at(2).mu);
  mu_errors.push_back(
      calibrated_cu_shield_signal_10_percent_20260114.peaks.at(2).mu_error);
  reduced_chi2.push_back(
      calibrated_cu_shield_signal_10_percent_20260114.reduced_chi2);

  //  run_names.push_back("Cu Shield Signal 90% (01/14/2026)");
  //  mu.push_back(calibrated_cu_shield_signal_90_percent_20260114.peaks.at(2).mu);
  //  mu_errors.push_back(
  //      calibrated_cu_shield_signal_90_percent_20260114.peaks.at(2).mu_error);
  //  reduced_chi2.push_back(
  //      calibrated_cu_shield_signal_90_percent_20260114.reduced_chi2);
  //
  //  run_names.push_back("No Shield Signal 5% (01/15/2026)");
  //  mu.push_back(calibrated_no_shield_signal_5_percent_20260115.peaks.at(1).mu);
  //  mu_errors.push_back(
  //      calibrated_no_shield_signal_5_percent_20260115.peaks.at(1).mu_error);
  //  reduced_chi2.push_back(
  //      calibrated_no_shield_signal_5_percent_20260115.reduced_chi2);
  //
  //  run_names.push_back("No Shield Graphite Castle Signal 10% (01/16/2026)");
  //  mu.push_back(
  //      calibrated_no_shield_graphite_castle_signal_10_percent_20260116.peaks.at(1).mu);
  //  mu_errors.push_back(
  //      calibrated_no_shield_graphite_castle_signal_10_percent_20260116.peaks.at(1)
  //          .mu_error);
  //  reduced_chi2.push_back(
  //      calibrated_no_shield_graphite_castle_signal_10_percent_20260116
  //          .reduced_chi2);

  std::cout << "Individual Run Results (Ge Peak mu):" << std::endl;

  for (size_t i = 0; i < mu.size(); ++i) {
    std::cout << std::left << std::setw(45) << run_names[i] << ": "
              << std::fixed << std::setprecision(4) << mu[i] << " +/- "
              << mu_errors[i] << " keV"
              << " (χ²/ndf = " << std::setprecision(3) << reduced_chi2[i] << ")"
              << std::endl;
  }

  Float_t sum_weights = 0.0;
  Float_t weighted_sum = 0.0;
  for (size_t i = 0; i < mu.size(); ++i) {
    Float_t weight = 1.0 / (mu_errors[i] * mu_errors[i]);
    weighted_sum += mu[i] * weight;
    sum_weights += weight;
  }
  Float_t combined_mu = weighted_sum / sum_weights;
  Float_t combined_error = std::sqrt(1.0 / sum_weights);

  std::cout << "Combined Ge mu: " << std::fixed << std::setprecision(4)
            << combined_mu << " +/- " << combined_error << " keV" << std::endl;
}
