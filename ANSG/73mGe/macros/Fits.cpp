#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "HyperEMGFitHelpers.hpp"
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
const Float_t E_BA133_81 = 80.9979;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;
const Float_t E_CD114M = 95.9023;

TH1F *LoadHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("calibrated_zoomedHist"));
  if (!hist) {
    std::cerr << "ERROR: Cannot find calibrated_zoomedHist in " << input_name
              << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }
  hist->SetDirectory(0);
  hist->Rebin(Constants::REBIN_FACTOR);
  hist->SetTitle(Form("; Energy [keV]; Counts / %d eV",
                      Constants::REBIN_FACTOR * Constants::BIN_WIDTH_EV));
  file->Close();
  delete file;
  return hist;
}

HyperEMGFitResult FitBackgroundPeak(const TString input_name,
                                    const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115) {
    fit_low = 67;
    fit_high = 77;
  } else if (input_name ==
             Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116) {
    fit_low = 67;
    fit_high = 80;
  }

  HyperEMGFitResult result = FitSingleHyperEMG(
      hist, fit_low, fit_high, kTRUE, input_name, "Background", interactive);
  delete hist;
  return result;
}

HyperEMGFitResult FitPbKAlpha(const TString input_name,
                              const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113) {
    fit_low = 66;
    fit_high = 81;
  } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114) {
    fit_low = 66;
    fit_high = 82;
  } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260113) {
    fit_low = 65;
    fit_high = 82;
  } else {
    fit_low = 65;
    fit_high = 81;
  }

  HyperEMGFitResult result =
      FitDoubleHyperEMG(hist, fit_low, fit_high, E_PB_KA1, E_PB_KA2, kTRUE,
                        input_name, "Pb_KAlpha", interactive);
  delete hist;
  return result;
}

HyperEMGFitResult FitSignalDoublePeak(const TString input_name,
                                      const HyperEMGPeakFitResult &bkg_peak,
                                      const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (input_name == Constants::NOSHIELDSIGNAL_5PERCENT_20260115) {
    fit_low = 64;
    fit_high = 77;
  } else if (input_name ==
             Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116) {
    fit_low = 60;
    fit_high = 77;
  } else {
    fit_low = 64;
    fit_high = 77;
  }

  // Double peak: Ge-73m (free) + constrained background peak
  Int_t ppp = FittingFunctions::kHyperEMGParamsPerPeak;
  Int_t npar = 2 * ppp + 2;
  TF1 *func =
      new TF1("hemg_Ge_" + input_name, &FittingFunctions::DoublePeakHyperEMG,
              fit_low, fit_high, npar);

  Double_t peak_height = hist->GetBinContent(hist->GetMaximumBin());

  // Peak 1: Ge-73m (free)
  FittingFunctions::SetHyperEMGParNames(func, 0, "Ge_");
  FittingFunctions::SetHyperEMGParLimits(func, 0);
  func->SetParLimits(0, fit_low, fit_high);
  FittingFunctions::SetHyperEMGInitialParams(func, 0, 68.75, 1.0,
                                             peak_height * 0.5);

  // Peak 2: constrained from background fit
  FittingFunctions::SetHyperEMGParNames(func, ppp, "Bkg_");
  func->FixParameter(ppp + 0, bkg_peak.mu);
  func->FixParameter(ppp + 1, bkg_peak.sigma);
  func->SetParameter(ppp + 2, bkg_peak.amplitude);
  func->SetParLimits(ppp + 2, 0, peak_height * 5);
  func->FixParameter(ppp + 3, bkg_peak.w_low1);
  func->FixParameter(ppp + 4, bkg_peak.tau_low1);
  func->FixParameter(ppp + 5, bkg_peak.w_low2);
  func->FixParameter(ppp + 6, bkg_peak.tau_low2);
  func->FixParameter(ppp + 7, bkg_peak.w_high1);
  func->FixParameter(ppp + 8, bkg_peak.tau_high1);
  func->FixParameter(ppp + 9, bkg_peak.w_high2);
  func->FixParameter(ppp + 10, bkg_peak.tau_high2);

  Int_t bkg_off = 2 * ppp;
  SetupHyperEMGBackground(func, bkg_off, peak_height, kTRUE);

  Double_t rlo = fit_low, rhi = fit_high;
  Bool_t fit_valid = kFALSE;
  Double_t chi2 = 0;

  HyperEMGFitResult result;
  result.peaks.emplace_back();
  result.peaks.emplace_back();

  if (interactive) {
    if (LoadHyperEMGParams(func, rlo, rhi, input_name, "Ge_double")) {
      TFitResultPtr refit = hist->Fit(func, "LSMRBENR+");
      if (refit.Get() && refit->IsValid()) {
        chi2 = refit->Chi2() / refit->Ndf();
        fit_valid = kTRUE;
      }
      std::cout << "Refit from saved params chi2/ndf = " << chi2 << std::endl;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(hist, func, rlo, rhi, 2,
                                     "Ge_double / " + input_name)) {
        func->GetRange(rlo, rhi);
        chi2 = func->GetChisquare() / func->GetNDF();
        SaveHyperEMGParams(func, rlo, rhi, input_name, "Ge_double");
        fit_valid = kTRUE;
        std::cout << "Interactive chi2/ndf = " << chi2 << std::endl;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    TFitResultPtr fit = hist->Fit(func, "LSMBNR");
    if (fit.Get() && fit->IsValid()) {
      chi2 = fit->Chi2() / fit->Ndf();
      fit_valid = kTRUE;
    }
  }

  if (fit_valid) {
    result.peaks[0] = ExtractHyperEMGPeak(func, 0);
    result.peaks[1] = ExtractHyperEMGPeak(func, ppp);
    result.bkg_constant = func->GetParameter(bkg_off);
    result.bkg_constant_error = func->GetParError(bkg_off);
    result.lin_bkg_slope = func->GetParameter(bkg_off + 1);
    result.lin_bkg_slope_error = func->GetParError(bkg_off + 1);
    result.reduced_chi2 = chi2;
    result.valid = kTRUE;
  }

  delete func;
  delete hist;
  return result;
}

HyperEMGFitResult FitSignalTriplePeak(const TString input_name,
                                      const HyperEMGFitResult &pb_fit,
                                      const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Double_t fit_low = 0, fit_high = 0;

  if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_20260113) {
    fit_low = 65;
    fit_high = 81;
  } else if (input_name == Constants::CDSHIELDSIGNAL_10PERCENT_20260113) {
    fit_low = 64;
    fit_high = 80;
  } else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_20260114) {
    fit_low = 62;
    fit_high = 80;
  } else if (input_name == Constants::CUSHIELDSIGNAL_90PERCENT_20260114) {
    fit_low = 63;
    fit_high = 80;
  } else {
    fit_low = 65;
    fit_high = 81;
  }

  // Triple peak: Ge-73m (free) + Pb-Ka1 + Pb-Ka2 (constrained from bkg fit)
  Int_t ppp = FittingFunctions::kHyperEMGParamsPerPeak;
  Int_t npar = 3 * ppp + 2;
  TF1 *func =
      new TF1("hemg_Ge_" + input_name, &FittingFunctions::TriplePeakHyperEMG,
              fit_low, fit_high, npar);

  Double_t peak_height = hist->GetBinContent(hist->GetMaximumBin());

  // Peak 1: Ge-73m (free)
  FittingFunctions::SetHyperEMGParNames(func, 0, "Ge_");
  FittingFunctions::SetHyperEMGParLimits(func, 0);
  func->SetParLimits(0, fit_low, fit_high);
  FittingFunctions::SetHyperEMGInitialParams(func, 0, 68.75, 1.0,
                                             peak_height * 0.3);

  // Peak 2: Pb-Ka1 (constrained shape, free amplitude)
  FittingFunctions::SetHyperEMGParNames(func, ppp, "PbKa1_");
  func->FixParameter(ppp + 0, pb_fit.peaks[0].mu);
  func->FixParameter(ppp + 1, pb_fit.peaks[0].sigma);
  func->SetParameter(ppp + 2, pb_fit.peaks[0].amplitude);
  func->SetParLimits(ppp + 2, 0, peak_height * 5);
  func->FixParameter(ppp + 3, pb_fit.peaks[0].w_low1);
  func->FixParameter(ppp + 4, pb_fit.peaks[0].tau_low1);
  func->FixParameter(ppp + 5, pb_fit.peaks[0].w_low2);
  func->FixParameter(ppp + 6, pb_fit.peaks[0].tau_low2);
  func->FixParameter(ppp + 7, pb_fit.peaks[0].w_high1);
  func->FixParameter(ppp + 8, pb_fit.peaks[0].tau_high1);
  func->FixParameter(ppp + 9, pb_fit.peaks[0].w_high2);
  func->FixParameter(ppp + 10, pb_fit.peaks[0].tau_high2);

  // Peak 3: Pb-Ka2 (constrained shape, free amplitude)
  FittingFunctions::SetHyperEMGParNames(func, 2 * ppp, "PbKa2_");
  func->FixParameter(2 * ppp + 0, pb_fit.peaks[1].mu);
  func->FixParameter(2 * ppp + 1, pb_fit.peaks[1].sigma);
  func->SetParameter(2 * ppp + 2, pb_fit.peaks[1].amplitude);
  func->SetParLimits(2 * ppp + 2, 0, peak_height * 5);
  func->FixParameter(2 * ppp + 3, pb_fit.peaks[1].w_low1);
  func->FixParameter(2 * ppp + 4, pb_fit.peaks[1].tau_low1);
  func->FixParameter(2 * ppp + 5, pb_fit.peaks[1].w_low2);
  func->FixParameter(2 * ppp + 6, pb_fit.peaks[1].tau_low2);
  func->FixParameter(2 * ppp + 7, pb_fit.peaks[1].w_high1);
  func->FixParameter(2 * ppp + 8, pb_fit.peaks[1].tau_high1);
  func->FixParameter(2 * ppp + 9, pb_fit.peaks[1].w_high2);
  func->FixParameter(2 * ppp + 10, pb_fit.peaks[1].tau_high2);

  Int_t bkg_off = 3 * ppp;
  SetupHyperEMGBackground(func, bkg_off, peak_height, kTRUE);

  Double_t rlo = fit_low, rhi = fit_high;
  Bool_t fit_valid = kFALSE;
  Double_t chi2 = 0;

  HyperEMGFitResult result;
  result.peaks.emplace_back();
  result.peaks.emplace_back();
  result.peaks.emplace_back();

  if (interactive) {
    if (LoadHyperEMGParams(func, rlo, rhi, input_name, "Ge_triple")) {
      TFitResultPtr refit = hist->Fit(func, "LSMRBENR+");
      if (refit.Get() && refit->IsValid()) {
        chi2 = refit->Chi2() / refit->Ndf();
        fit_valid = kTRUE;
      }
      std::cout << "Refit from saved params chi2/ndf = " << chi2 << std::endl;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(hist, func, rlo, rhi, 3,
                                     "Ge_triple / " + input_name)) {
        func->GetRange(rlo, rhi);
        chi2 = func->GetChisquare() / func->GetNDF();
        SaveHyperEMGParams(func, rlo, rhi, input_name, "Ge_triple");
        fit_valid = kTRUE;
        std::cout << "Interactive chi2/ndf = " << chi2 << std::endl;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    TFitResultPtr fit = hist->Fit(func, "LSMBNR");
    if (fit.Get() && fit->IsValid()) {
      chi2 = fit->Chi2() / fit->Ndf();
      fit_valid = kTRUE;
    }
  }

  if (fit_valid) {
    result.peaks[0] = ExtractHyperEMGPeak(func, 0);
    result.peaks[1] = ExtractHyperEMGPeak(func, ppp);
    result.peaks[2] = ExtractHyperEMGPeak(func, 2 * ppp);
    result.bkg_constant = func->GetParameter(bkg_off);
    result.bkg_constant_error = func->GetParError(bkg_off);
    result.lin_bkg_slope = func->GetParameter(bkg_off + 1);
    result.lin_bkg_slope_error = func->GetParError(bkg_off + 1);
    result.reduced_chi2 = chi2;
    result.valid = kTRUE;
  }

  delete func;
  delete hist;
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

  HyperEMGFitResult cd_bkg_10 = FitPbKAlpha(
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113, interactive);
  HyperEMGFitResult cd_bkg_25 = FitPbKAlpha(
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113, interactive);
  HyperEMGFitResult cu_bkg_0113 = FitPbKAlpha(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113, interactive);
  HyperEMGFitResult cu_bkg_0114 = FitPbKAlpha(
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114, interactive);

  // Triple peak signal fits (Ge peak + constrained Pb K-alpha)

  HyperEMGFitResult cd_sig_10 = FitSignalTriplePeak(
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113, cd_bkg_10, interactive);
  HyperEMGFitResult cd_sig_25 = FitSignalTriplePeak(
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113, cd_bkg_25, interactive);
  HyperEMGFitResult cu_sig_0113 = FitSignalTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113, cu_bkg_0113, interactive);
  HyperEMGFitResult cu_sig_0114 = FitSignalTriplePeak(
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114, cu_bkg_0114, interactive);

  if (cd_sig_10.valid) {
    run_names.push_back("Cd Shield Signal 10% (01/13)");
    mu.push_back(cd_sig_10.peaks.at(0).mu);
    mu_errors.push_back(cd_sig_10.peaks.at(0).mu_error);
    reduced_chi2.push_back(cd_sig_10.reduced_chi2);
  }

  if (cd_sig_25.valid) {
    run_names.push_back("Cd Shield Signal 25% (01/13)");
    mu.push_back(cd_sig_25.peaks.at(0).mu);
    mu_errors.push_back(cd_sig_25.peaks.at(0).mu_error);
    reduced_chi2.push_back(cd_sig_25.reduced_chi2);
  }

  if (cu_sig_0113.valid) {
    run_names.push_back("Cu Shield Signal 10% (01/13)");
    mu.push_back(cu_sig_0113.peaks.at(0).mu);
    mu_errors.push_back(cu_sig_0113.peaks.at(0).mu_error);
    reduced_chi2.push_back(cu_sig_0113.reduced_chi2);
  }

  if (cu_sig_0114.valid) {
    run_names.push_back("Cu Shield Signal 10% (01/14)");
    mu.push_back(cu_sig_0114.peaks.at(0).mu);
    mu_errors.push_back(cu_sig_0114.peaks.at(0).mu_error);
    reduced_chi2.push_back(cu_sig_0114.reduced_chi2);
  }

  // Double peak signal fits (Ge + constrained bkg peak, no Pb overlap)

  HyperEMGFitResult noshield_bkg = FitBackgroundPeak(
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115, interactive);
  if (noshield_bkg.valid) {
    HyperEMGFitResult noshield_sig =
        FitSignalDoublePeak(Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                            noshield_bkg.peaks[0], interactive);
    if (noshield_sig.valid) {
      run_names.push_back("No Shield Signal 5% (01/15)");
      mu.push_back(noshield_sig.peaks.at(0).mu);
      mu_errors.push_back(noshield_sig.peaks.at(0).mu_error);
      reduced_chi2.push_back(noshield_sig.reduced_chi2);
    }
  }

  HyperEMGFitResult graphite_bkg = FitBackgroundPeak(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      interactive);
  if (graphite_bkg.valid) {
    HyperEMGFitResult graphite_sig = FitSignalDoublePeak(
        Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
        graphite_bkg.peaks[0], interactive);
    if (graphite_sig.valid) {
      run_names.push_back("No Shield Graphite Castle 10% (01/16)");
      mu.push_back(graphite_sig.peaks.at(0).mu);
      mu_errors.push_back(graphite_sig.peaks.at(0).mu_error);
      reduced_chi2.push_back(graphite_sig.reduced_chi2);
    }
  }

  // Summary

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
