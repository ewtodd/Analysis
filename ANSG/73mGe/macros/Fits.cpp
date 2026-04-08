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

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_81 = 80.98;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;
const Float_t E_CD114M = 95.9023;

TH1D *LoadHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1D *hist = static_cast<TH1D *>(file->Get("calibrated_zoomedHist"));
  if (!hist) {
    std::cerr << "ERROR: Cannot find calibrated_zoomedHist in " << input_name
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

FitResult FitCalibrationPeak(const TString input_name, const TString peak_name,
                             const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (peak_name == "Am_59.5keV") {
    if (input_name == Constants::POSTREACTOR_AM241_20260113)
      fitter = new FittingUtils(hist, 50, 70, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else if (input_name == Constants::POSTREACTOR_AM241_BA133_20260116)
      fitter = new FittingUtils(hist, 55, 70, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else
      fitter = new FittingUtils(hist, 51, 71, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
  }
  if (peak_name == "Ba_80.98keV") {
    fitter =
        new FittingUtils(hist, 75, 90, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  }
  if (peak_name == "Cd114m_95.9keV") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115)
      fitter = new FittingUtils(hist, 91, 100, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
    else
      fitter = new FittingUtils(hist, 91, 103, use_flat_background, use_step,
                                use_low_exp_tail, use_low_lin_tail,
                                use_high_exp_tail);
  }

  if (interactive)
    fitter->SetInteractive();
  FitResult result = fitter->FitSinglePeak(input_name, peak_name);
  delete hist;
  delete fitter;
  return result;
}

FitResult FitBackgroundPeak(const TString input_name,
                            const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115)
    fitter =
        new FittingUtils(hist, 67, 77, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name ==
           Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116)
    fitter =
        new FittingUtils(hist, 67, 80, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter->SetInteractive();
  FitResult result = fitter->FitSinglePeak(input_name, "Background");
  delete hist;
  delete fitter;
  return result;
}

FitResult FitPbKAlpha(const TString input_name, const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113)
    fitter =
        new FittingUtils(hist, 66, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114)
    fitter =
        new FittingUtils(hist, 66, 82, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260113) {
    use_flat_background = kFALSE;
    fitter =
        new FittingUtils(hist, 65, 82, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  } else {
    use_flat_background = kFALSE;
    fitter =
        new FittingUtils(hist, 65, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  }

  if (interactive)
    fitter->SetInteractive();
  FitResult result =
      fitter->FitDoublePeak(input_name, "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
  delete hist;
  delete fitter;
  return result;
}

FitResult FitSignalDoublePeak(const TString input_name,
                              const PeakFitResult &constrained_peak,
                              const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (input_name == Constants::NOSHIELDSIGNAL_5PERCENT_20260115)
    fitter =
        new FittingUtils(hist, 64, 77, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name ==
           Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116)
    fitter =
        new FittingUtils(hist, 60, 77, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter->SetInteractive();
  FitResult result =
      fitter->FitDoublePeak(input_name, "Ge", constrained_peak, 68.75);
  delete hist;
  delete fitter;
  return result;
}

FitResult FitSignalTriplePeak(const TString input_name,
                              const FitResult &constrained_peaks,
                              const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_20260113)
    fitter =
        new FittingUtils(hist, 65, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CDSHIELDSIGNAL_10PERCENT_20260113)
    fitter =
        new FittingUtils(hist, 64, 80, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_20260114)
    fitter =
        new FittingUtils(hist, 62, 80, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else if (input_name == Constants::CUSHIELDSIGNAL_90PERCENT_20260114)
    fitter =
        new FittingUtils(hist, 63, 80, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  else
    fitter =
        new FittingUtils(hist, 65, 81, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter->SetInteractive();
  FitResult result =
      fitter->FitTriplePeak(input_name, "Ge", constrained_peaks, 68.75);
  delete hist;
  delete fitter;
  return result;
}

void Fits() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  std::vector<TString> run_names;
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> reduced_chi2;

  // Pb K-alpha backgrounds (used as constraints for signal fits)

  FitResult cd_bkg_10 = FitPbKAlpha(
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113, interactive);
  FitResult cd_bkg_25 = FitPbKAlpha(
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113, interactive);
  FitResult cu_bkg_0113 = FitPbKAlpha(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113, interactive);
  FitResult cu_bkg_0114 = FitPbKAlpha(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, interactive);

  // Triple peak signal fits (Ge peak + constrained Pb K-alpha)

  FitResult cd_sig_10 = FitSignalTriplePeak(
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113, cd_bkg_10, interactive);
  FitResult cd_sig_25 = FitSignalTriplePeak(
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113, cd_bkg_25, interactive);
  FitResult cu_sig_0113 = FitSignalTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113, cu_bkg_0113, interactive);
  FitResult cu_sig_0114 = FitSignalTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, cu_bkg_0114, interactive);

  run_names.push_back("Cd Shield Signal 10% (01/13)");
  mu.push_back(cd_sig_10.peaks.at(0).mu);
  mu_errors.push_back(cd_sig_10.peaks.at(0).mu_error);
  reduced_chi2.push_back(cd_sig_10.reduced_chi2);

  run_names.push_back("Cd Shield Signal 25% (01/13)");
  mu.push_back(cd_sig_25.peaks.at(0).mu);
  mu_errors.push_back(cd_sig_25.peaks.at(0).mu_error);
  reduced_chi2.push_back(cd_sig_25.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% (01/13)");
  mu.push_back(cu_sig_0113.peaks.at(0).mu);
  mu_errors.push_back(cu_sig_0113.peaks.at(0).mu_error);
  reduced_chi2.push_back(cu_sig_0113.reduced_chi2);

  run_names.push_back("Cu Shield Signal 10% (01/14)");
  mu.push_back(cu_sig_0114.peaks.at(0).mu);
  mu_errors.push_back(cu_sig_0114.peaks.at(0).mu_error);
  reduced_chi2.push_back(cu_sig_0114.reduced_chi2);

  std::cout << std::endl;
  std::cout << "Individual Run Results (Ge Peak mu):" << std::endl;

  for (size_t i = 0; i < mu.size(); ++i) {
    std::cout << std::left << std::setw(50) << run_names[i] << ": "
              << std::fixed << std::setprecision(4) << mu[i] << " +/- "
              << mu_errors[i] << " keV"
              << " (chi2/ndf = " << std::setprecision(3) << reduced_chi2[i]
              << ")" << std::endl;
  }

  Float_t sum_weights = 0.0;
  Float_t weighted_sum = 0.0;
  for (size_t i = 0; i < mu.size(); ++i) {
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
    std::cout << "Combined Ge mu: " << std::fixed << std::setprecision(4)
              << combined_mu << " +/- " << combined_error << " keV"
              << std::endl;
  }
}
