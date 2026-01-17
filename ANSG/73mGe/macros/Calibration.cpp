#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <vector>

FitResultGaussianTailStep
FitSinglePeak(const TString input_name, const TString peak_name,
              const TString branch_name, const TString tree_name,
              const TString formatted_branch_name, const Int_t color,
              const Float_t expected_mu) {

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

  FittingUtils *fitter = new FittingUtils(kFALSE, kTRUE, kFALSE);
  fitter->LoadProcessed(input_name, branch_name, tree_name);
  fitter->SetNumHistBins(Constants::ZOOMED_NBINS);
  fitter->SetMinHistValue(Constants::ZOOMED_XMIN);
  fitter->SetMaxHistValue(Constants::ZOOMED_XMAX);
  fitter->SetExpectedMu(expected_mu);
  fitter->SetFitRange(expected_mu - 0.20 * expected_mu,
                      expected_mu + 0.12 * expected_mu);
  FitResultGaussianTailStep result;

  if (peak_name == "Am_59.5keV") {
    fitter->SetExpectedAmplitude(11000);
    fitter->SetExpectedTail(10);
    fitter->SetExpectedTailAmplitude(200);
    fitter->SetExpectedStepAmplitude(200);
    fitter->SetExpectedMu(59.6788);
    fitter->SetExpectedSigma(1);
    fitter->GetFitFunction()->SetParLimits(0, 0, 20000);
  }

  result = fitter->FitPeakGaussianLowTailLowStep(canvas, color, peak_name,
                                                 formatted_branch_name);
  return result;
  delete fitter;
  delete canvas;
}

std::vector<FitResultGaussianTailStep> FitMultiplePeaks(
    std::vector<TString> input_names, std::vector<TString> peak_names,
    std::vector<Float_t> mu_guesses, const TString formatted_branch_name,
    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {

  std::vector<FitResultGaussianTailStep> results;
  Int_t entries = input_names.size();
  for (Int_t i = 0; i < entries; i++) {
    FitResultGaussianTailStep result;
    if (input_names[i] != "zero") {

      Bool_t isFiltered = input_names[i].Contains("_filtered");

      if (isFiltered && !Constants::FILTERED) {
        std::cerr << "File name contains filtered but constant is not set."
                  << std::endl;
      }

      TString tree_name;
      TString energyBranchName;

      if (isFiltered) {
        tree_name = "bef_tree";
        energyBranchName = "energy";
      } else {
        tree_name = "bef_tree_event_summary";
        energyBranchName = "totalEnergy";
      }

      result = FitSinglePeak(input_names[i], peak_names[i], energyBranchName,
                             tree_name, formatted_branch_name, colors[i],
                             mu_guesses[i]);
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

TF1 *CreateAndSaveCalibration(std::vector<Float_t> mu,
                              std::vector<Float_t> calibration_values_keV,
                              std::vector<Float_t> mu_errors) {

  Int_t size = calibration_values_keV.size();

  TGraphErrors *calibration_curve =
      new TGraphErrors(size, mu.data(), calibration_values_keV.data(),
                       mu_errors.data(), nullptr);
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(
      calibration_curve, kBlue,
      "; Precalibrated Energy [keV]; Deposited Energy [keV]");

  Double_t x_min, x_max, y_min, y_max;
  calibration_curve->ComputeRange(x_min, y_min, x_max, y_max);
  calibration_curve->GetXaxis()->SetRangeUser(0, x_max * 1.2);
  calibration_curve->GetYaxis()->SetRangeUser(0, y_max * 1.2);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(21);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("APE");

  TF1 *calibration_fit = new TF1("linear", "pol1", 0, 5000);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParLimits(0, -1e-2, 1e-2);
  calibration_fit->SetParameter(1, 0.076);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit);
  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);
  return calibration_fit;
}

void PulseHeightToDepositedEnergy(
    std::vector<TString> input_names, TF1 *calibration_function,
    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors()) {
  Int_t entries = input_names.size();

  TString calibration_function_filepath =
      "root_files/calibration_function.root";
  TFile *calibration_file =
      new TFile(calibration_function_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);

  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];
    TH1F *light_output_hist =
        new TH1F("", "; Light Output [keVee]; Counts", 500, 0,
                 calibration_function->Eval(16384));
    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");

    output->cd();

    TTree *features_tree = static_cast<TTree *>(output->Get("features"));

    Float_t pulse_height, light_output_keVee;
    features_tree->SetBranchAddress("pulse_height", &pulse_height);
    features_tree->Branch("light_output", &light_output_keVee,
                          "light_output/F");

    Int_t num_entries = features_tree->GetEntries();

    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      light_output_keVee = calibration_function->Eval(pulse_height);
      features_tree->GetBranch("light_output")->Fill();
      light_output_hist->Fill(light_output_keVee);
    }

    PlottingUtils::ConfigureAndDrawHistogram(light_output_hist, colors[i]);
    PlottingUtils::SaveFigure(canvas, input_name + "_light_output.png");

    features_tree->Write("", TObject::kOverwrite);
    light_output_hist->Write("Light Output", TObject::kOverwrite);
    output->Close();
    delete canvas;
    delete light_output_hist;
  }
  calibration_file->Close();
}

void Calibration() {
  InitUtils::SetROOTPreferences();

  std::vector<Float_t> calibration_values_keV = {0,
                                                 59.5409, // 72.8042, 74.9694,
                                                 84.936};
  std::vector<Float_t> mu_guesses = {0, 59.5409, // 72.8042, 74.9694,
                                     84.936};

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  TString input_name_Am241 = "01132026-PostReactor-Calibration" + suffix;
  TString input_name_PbXRays =
      "01142026-CuShield-ActiveBackground-10Percent" + suffix;

  std::vector<TString> input_names_calibrations = {"zero", input_name_Am241,
                                                   input_name_PbXRays};

  std::vector<TString> peak_names = {"zero",
                                     "Am_59.5keV", //"Pb_72.8keV", "Pb_75.0keV",
                                     "Pb_84.9keV"};

  std::vector<Int_t> defaultColors = PlottingUtils::GetDefaultColors();
  std::vector<Int_t> colors = {0,
                               defaultColors[1],
                               defaultColors[0],
                               defaultColors[1],
                               defaultColors[1],
                               defaultColors[1]};

  std::vector<FitResultGaussianTailStep> fit_results =
      FitMultiplePeaks(input_names_calibrations, peak_names, mu_guesses,
                       "Precalibrated Energy [keV]", colors);

  Int_t size = calibration_values_keV.size();

  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;

  for (Int_t i = 0; i < size; i++) {
    mu.push_back(fit_results[i].mu);
    mu_errors.push_back(fit_results[i].mu_error);
  }

  //  TF1 *calibration_function =
  //      CreateAndSaveCalibration(mu, calibration_values_keV, mu_errors);
  //
  //  std::vector<TString> input_names = {input_name_Am241, input_name_PbXRays};
  //  PulseHeightToDepositedEnergy(input_names, calibration_function);
}
