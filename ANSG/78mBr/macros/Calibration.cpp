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

struct CalibrationData {
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> peak_names;
};

FitResultStandard FitSinglePeak(const TString input_name,
                                const TString peak_name,
                                const Float_t expected_mu) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << input_name << ".root" << std::endl;
    return {};
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("Pulse Height"));
  if (!hist) {
    std::cerr << "Cannot find 'Pulse Height' histogram in " << input_name
              << ".root" << std::endl;
    file->Close();
    delete file;
    return {};
  }

  hist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultStandard result = {};

  Float_t fit_low = expected_mu * 0.85;
  Float_t fit_high = expected_mu * 1.15;

  if (peak_name == "La_33keV") {
    fit_low = 250;
    fit_high = 800;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Am_59keV") {
    fit_low = 600;
    fit_high = 900;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Eu_122keV") {
    fit_low = expected_mu * 0.90;
    fit_high = expected_mu * 1.08;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Eu_245keV") {
    fit_low = 2990;
    fit_high = 3500;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Eu_344keV") {
    fit_low = 4100;
    fit_high = 5000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else {
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  }

  result = fitter->FitPeakStandard(input_name, peak_name);

  delete hist;
  delete fitter;
  return result;
}

CalibrationData FitCalibrationPeaks() {
  CalibrationData cal_data;

  // Zero point
  cal_data.peak_names.push_back("Zero");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);

  // La K-alpha 33 keV
  FitResultStandard la_result = FitSinglePeak(
      Constants::CALIBRATION_EU152, "La_33keV", 525); // ~525 ADC expected
  cal_data.peak_names.push_back("La_33keV");
  cal_data.mu.push_back(la_result.mu);
  cal_data.mu_errors.push_back(la_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_LA_33KEV);
  cal_data.reduced_chi2.push_back(la_result.reduced_chi2);

  // Am-241 59.5 keV
  FitResultStandard am_result = FitSinglePeak(
      Constants::CALIBRATION_AM241, "Am_59keV", 780); // ~780 ADC expected
  cal_data.peak_names.push_back("Am_59keV");
  cal_data.mu.push_back(am_result.mu);
  cal_data.mu_errors.push_back(am_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_AM241_59KEV);
  cal_data.reduced_chi2.push_back(am_result.reduced_chi2);

  // Eu-152 122 keV
  FitResultStandard eu122_result = FitSinglePeak(
      Constants::CALIBRATION_EU152, "Eu_122keV", 1650); // ~1650 ADC expected
  cal_data.peak_names.push_back("Eu_122keV");
  cal_data.mu.push_back(eu122_result.mu);
  cal_data.mu_errors.push_back(eu122_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_EU152_122KEV);
  cal_data.reduced_chi2.push_back(eu122_result.reduced_chi2);

  // Eu-152 245 keV
  FitResultStandard eu245_result = FitSinglePeak(
      Constants::CALIBRATION_EU152, "Eu_245keV", 3275); // ~3250 ADC expected
  cal_data.peak_names.push_back("Eu_245keV");
  cal_data.mu.push_back(eu245_result.mu);
  cal_data.mu_errors.push_back(eu245_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_EU152_245KEV);
  cal_data.reduced_chi2.push_back(eu245_result.reduced_chi2);

  // Eu-152 344 keV
  FitResultStandard eu344_result = FitSinglePeak(
      Constants::CALIBRATION_EU152, "Eu_344keV", 4577); // ~4577 ADC expected
  cal_data.peak_names.push_back("Eu_344keV");
  cal_data.mu.push_back(eu344_result.mu);
  cal_data.mu_errors.push_back(eu344_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_EU152_344KEV);
  cal_data.reduced_chi2.push_back(eu344_result.reduced_chi2);

  return cal_data;
}

void PrintCalibrationSummary(const CalibrationData &cal_data) {
  for (size_t i = 0; i < cal_data.mu.size(); ++i) {
    std::cout << cal_data.peak_names[i] << ": " << std::fixed
              << std::setprecision(2) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " ADC -> "
              << cal_data.calibration_values_keV[i] << " keV" << std::endl;
  }
}

TF1 *CreateAndSaveCalibration(const CalibrationData &cal_data) {

  Int_t size = cal_data.calibration_values_keV.size();

  TGraphErrors *calibration_curve = new TGraphErrors(
      size, cal_data.mu.data(), cal_data.calibration_values_keV.data(),
      cal_data.mu_errors.data(), nullptr);

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
  calibration_fit->SetParLimits(0, -1e-2, 1e-2);
  calibration_fit->SetParameter(1, 0.076);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");
  calibration_fit->Draw("SAME");
  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);

  delete canvas;
  return calibration_fit;
}

void PulseHeightToLightOutput(const std::vector<TString> &input_names,
                              TF1 *calibration_function) {

  TString calibration_function_filepath =
      "root_files/calibration_function.root";
  TFile *calibration_file =
      new TFile(calibration_function_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);
  calibration_file->Close();
  delete calibration_file;

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  for (size_t i = 0; i < input_names.size(); i++) {
    TString input_name = input_names[i];
    Int_t color = colors[i % colors.size()];

    TH1F *light_output_hist =
        new TH1F("",
                 Form("; Light Output [keVee]; Counts / %.1d keV",
                      Constants::LO_BIN_WIDTH),
                 Constants::LO_HIST_NBINS, Constants::LO_HIST_XMIN,
                 Constants::LO_HIST_XMAX);

    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);

    TString output_filepath = "root_files/" + input_name + ".root";
    TFile *output = new TFile(output_filepath, "UPDATE");

    if (!output || output->IsZombie()) {
      std::cerr << "Cannot open " << output_filepath << std::endl;
      delete light_output_hist;
      delete canvas;
      continue;
    }

    output->cd();

    TTree *features_tree = static_cast<TTree *>(output->Get("features"));
    if (!features_tree) {
      std::cerr << "Cannot find 'features' tree in " << output_filepath
                << std::endl;
      output->Close();
      delete output;
      delete light_output_hist;
      delete canvas;
      continue;
    }

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

    PlottingUtils::ConfigureAndDrawHistogram(light_output_hist, color);
    PlottingUtils::SaveFigure(canvas, input_name + "_light_output.png");

    features_tree->Write("", TObject::kOverwrite);
    light_output_hist->Write("Light Output", TObject::kOverwrite);
    output->Close();
    delete output;
    delete canvas;
    delete light_output_hist;
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences();

  CalibrationData cal_data = FitCalibrationPeaks();

  PrintCalibrationSummary(cal_data);

  TF1 *calibration_function = CreateAndSaveCalibration(cal_data);

  std::vector<TString> datasets_to_calibrate = {
      Constants::CALIBRATION_AM241, Constants::CALIBRATION_EU152,
      Constants::BACKGROUND,        Constants::IRRADIATION_ONE,
      Constants::IRRADIATION_TWO,   Constants::IRRADIATION_THREE,
      Constants::IRRADIATION_FOUR};

  PulseHeightToLightOutput(datasets_to_calibrate, calibration_function);
}
