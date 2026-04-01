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

static const Int_t N_CRYSTALS = 4;

TH1F *LoadHistogram(const TString input_name, Int_t crystal,
                    const Bool_t use_calibrated = kFALSE) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TString histName = use_calibrated
      ? Form("calibrated_zoomedHist_crystal%d", crystal)
      : Form("zoomedHist_crystal%d", crystal);
  TH1F *hist = static_cast<TH1F *>(file->Get(histName));
  if (!hist) {
    std::cerr << "ERROR: Cannot find " << histName << " in " << input_name
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

// Calibration single peaks
FitResult FitCalibrationPeak(const TString input_name, const TString peak_name,
                             Int_t crystal, const Bool_t use_calibrated,
                             const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
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
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result = fitter->FitSinglePeak(label, peak_name);
  delete hist;
  delete fitter;
  return result;
}

// Background single peaks (for constrained signal fits)
FitResult FitBackgroundPeak(const TString input_name, Int_t crystal,
                            const Bool_t use_calibrated,
                            const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
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
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result = fitter->FitSinglePeak(label, "Background");
  delete hist;
  delete fitter;
  return result;
}

// Pb K-alpha double peak (shared between calibration and signal analysis)
FitResult FitPbKAlpha(const TString input_name, Int_t crystal,
                      const Bool_t use_calibrated,
                      const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
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
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result =
      fitter->FitDoublePeak(label, "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
  delete hist;
  delete fitter;
  return result;
}

// Constrained double peak (signal with one constrained peak from background)
FitResult FitSignalDoublePeak(const TString input_name,
                              const PeakFitResult &constrained_peak,
                              Int_t crystal, const Bool_t use_calibrated,
                              const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
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
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result =
      fitter->FitDoublePeak(label, "Ge", constrained_peak, 68.75);
  delete hist;
  delete fitter;
  return result;
}

// Triple peak (signal with two constrained peaks from background)
FitResult FitSignalTriplePeak(const TString input_name,
                              const FitResult &constrained_peaks,
                              Int_t crystal, const Bool_t use_calibrated,
                              const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, crystal, use_calibrated);
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
  TString label = Form("%s_crystal%d", input_name.Data(), crystal);
  FitResult result =
      fitter->FitTriplePeak(label, "Ge", constrained_peaks, 68.75);
  delete hist;
  delete fitter;
  return result;
}

void Fits() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t use_calibrated = kTRUE;
  Bool_t interactive = kTRUE;

  std::vector<TString> run_names;
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> reduced_chi2;

  for (Int_t c = 0; c < N_CRYSTALS; c++) {
    std::cout << std::endl << "=== Crystal " << c << " ===" << std::endl;

    // Pb K-alpha backgrounds (used as constraints for signal fits)

    FitResult cd_bkg_10 =
        FitPbKAlpha(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
                    c, use_calibrated, interactive);
    FitResult cd_bkg_25 =
        FitPbKAlpha(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
                    c, use_calibrated, interactive);
    FitResult cu_bkg_0113 =
        FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                    c, use_calibrated, interactive);
    FitResult cu_bkg_0114 =
        FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                    c, use_calibrated, interactive);

    // Triple peak signal fits (Ge peak + constrained Pb K-alpha)

    FitResult cd_sig_10 =
        FitSignalTriplePeak(Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                            cd_bkg_10, c, use_calibrated, interactive);
    FitResult cd_sig_25 =
        FitSignalTriplePeak(Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                            cd_bkg_25, c, use_calibrated, interactive);
    FitResult cu_sig_0113 =
        FitSignalTriplePeak(Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                            cu_bkg_0113, c, use_calibrated, interactive);
    FitResult cu_sig_0114 =
        FitSignalTriplePeak(Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                            cu_bkg_0114, c, use_calibrated, interactive);

    run_names.push_back(Form("Cd Shield Signal 10%% (01/13) crystal %d", c));
    mu.push_back(cd_sig_10.peaks.at(0).mu);
    mu_errors.push_back(cd_sig_10.peaks.at(0).mu_error);
    reduced_chi2.push_back(cd_sig_10.reduced_chi2);

    run_names.push_back(Form("Cd Shield Signal 25%% (01/13) crystal %d", c));
    mu.push_back(cd_sig_25.peaks.at(0).mu);
    mu_errors.push_back(cd_sig_25.peaks.at(0).mu_error);
    reduced_chi2.push_back(cd_sig_25.reduced_chi2);

    run_names.push_back(Form("Cu Shield Signal 10%% (01/13) crystal %d", c));
    mu.push_back(cu_sig_0113.peaks.at(0).mu);
    mu_errors.push_back(cu_sig_0113.peaks.at(0).mu_error);
    reduced_chi2.push_back(cu_sig_0113.reduced_chi2);

    run_names.push_back(Form("Cu Shield Signal 10%% (01/14) crystal %d", c));
    mu.push_back(cu_sig_0114.peaks.at(0).mu);
    mu_errors.push_back(cu_sig_0114.peaks.at(0).mu_error);
    reduced_chi2.push_back(cu_sig_0114.reduced_chi2);
  }

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
