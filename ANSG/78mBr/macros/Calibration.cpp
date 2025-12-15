#include "FitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

FitResult FitSinglePeak(const TString input_name, const TString peak_name,
                        const TString branch_name,
                        const TString formatted_branch_name, const Int_t color,
                        const Float_t expected_mu) {
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

  FitUtils *fitter = new FitUtils();
  fitter->LoadProcessed(input_name, branch_name);
  fitter->SetNumHistBins(1500);
  fitter->SetMaxHistValue(16384);
  fitter->SetExpectedMu(expected_mu);
  fitter->SetFitRange(expected_mu - 0.1 * expected_mu,
                      expected_mu + 0.1 * expected_mu);
  FitResult result;

  if (peak_name == "Ba_37keV") {
    fitter->SetExpectedAmplitude(3500);
    fitter->SetExpectedSigma(50);
    fitter->SetFitRange(expected_mu - 0.18 * expected_mu,
                        expected_mu + 0.18 * expected_mu);
  }

  if (peak_name == "Am_59keV") {
    fitter->SetExpectedAmplitude(550);
    fitter->SetExpectedSigma(34);
    fitter->SetFitRange(expected_mu - 0.13 * expected_mu,
                        expected_mu + 0.13 * expected_mu);
    TF1 *function = fitter->GetFitFunction();
    function->SetParLimits(0, 0, 3000);
  }

  if (peak_name == "Eu_122keV") {
    fitter->SetExpectedSigma(50);
    fitter->SetFitRange(expected_mu - 0.1 * expected_mu,
                        expected_mu + 0.08 * expected_mu);
  }

  if (peak_name == "Eu_344keV") {
    fitter->SetExpectedAmplitude(1600);
    fitter->SetExpectedSigma(187);
    fitter->SetFitRange(expected_mu - 0.12 * expected_mu,
                        expected_mu + 0.12 * expected_mu);
  }

  if (peak_name == "Eu_411keV") {
    fitter->SetExpectedSigma(100);
    TF1 *function = fitter->GetFitFunction();
    function->SetParLimits(0, 0, 3000);
    function->SetParLimits(2, 40, 300);
    fitter->SetNumHistBins(500);
    function->SetParLimits(4, 0, 1e2);
  }
  if (peak_name == "Eu_444eV") {
    fitter->SetExpectedSigma(250);
    TF1 *function = fitter->GetFitFunction();
    function->SetParLimits(0, 0, 3000);
    function->SetParLimits(2, 100, 300);
    fitter->SetNumHistBins(500);
    function->SetParLimits(4, -1e2, 0);
  }

  result = fitter->FitPeak(canvas, color, peak_name, formatted_branch_name);
  return result;
  delete fitter;
  delete canvas;
}

std::vector<FitResult> FitMultiplePeaks(
    std::vector<TString> input_names, std::vector<TString> peak_names,
    std::vector<Float_t> mu_guesses, const TString branch_name,
    const TString formatted_branch_name,
    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {

  std::vector<FitResult> results;
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    FitResult result;
    if (input_names[i] != "zero") {
      result = FitSinglePeak(input_names[i], peak_names[i], branch_name,
                             formatted_branch_name, colors[i], mu_guesses[i]);
    } else {
      result.mu = 0;
      result.mu_error = 0;
      result.sigma = 0;
      result.sigma_error = 0;
    }
    results.push_back(result);
  }
  return results;
}

void Calibration() {
  PlottingUtils::SetROOTPreferences();

  std::vector<Float_t> calibration_values_keV = {0,     37.4,  59.5409,
                                                 121.8, 244.7, 344.3};
  std::vector<Float_t> calibration_uncertainties_keV = {0, 0, 0, 0, 0, 0};
  std::vector<Float_t> mu_guesses = {0, 479, 780, 1650, 3250, 4577};

  TString input_name_Am241 = "calibration_Am241";
  TString input_name_Eu152 = "calibration_Eu152";
  TString input_name_bkg = "background";

  std::vector<TString> input_names = {"zero",           input_name_bkg,
                                      input_name_Am241, input_name_Eu152,
                                      input_name_Eu152, input_name_Eu152};

  std::vector<TString> peak_names = {"zero",      "Ba_37keV",  "Am_59keV",
                                     "Eu_122keV", "Eu_244keV", "Eu_344keV"};

  std::vector<Int_t> defaultColors = PlottingUtils::GetDefaultColors();
  std::vector<Int_t> colors = {0,
                               defaultColors[2],
                               defaultColors[0],
                               defaultColors[1],
                               defaultColors[1],
                               defaultColors[1]};

  std::vector<FitResult> fit_results =
      FitMultiplePeaks(input_names, peak_names, mu_guesses, "pulse_height",
                       "Pulse Height [ADC]", colors);

  Int_t size = calibration_values_keV.size();

  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;

  for (Int_t i = 0; i < size; i++) {
    mu.push_back(fit_results[i].mu);
    mu_errors.push_back(fit_results[i].mu_error);
  }

  TGraphErrors *calibration_curve =
      new TGraphErrors(size, mu.data(), calibration_values_keV.data(),
                       mu_errors.data(), calibration_uncertainties_keV.data());
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(calibration_curve, kBlue,
                                "; Pulse Height [ADC]; Deposited Energy [keV]");

  Double_t x_min, x_max, y_min, y_max;
  calibration_curve->ComputeRange(x_min, y_min, x_max, y_max);
  calibration_curve->GetXaxis()->SetRangeUser(-250, x_max * 1.2);
  calibration_curve->GetYaxis()->SetRangeUser(-25, y_max * 1.2);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(21);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("APE");

  TF1 *calibration_fit = new TF1("linear", "pol1", 0, 5000);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParLimits(0, -1e-1, 1e-1);
  calibration_fit->SetParameter(1, 0.076);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit);
  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);
}
