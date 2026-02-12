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

FitResultDetailed FitSinglePeakDetailed(const TString input_name,
                                        const TString peak_name,
                                        const Float_t expected_mu) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return {};
  } else
    std::cout << "FOUND " << input_name << ".root" << std::endl;

  TH1F *hist = static_cast<TH1F *>(file->Get("long_integral"));
  if (!hist) {
    std::cerr << "Cannot find 'Pulse Integral' histogram in " << input_name
              << ".root" << std::endl;
    file->Close();
    delete file;
    return {};
  }

  hist->SetDirectory(0);
  file->Close();
  delete file;

  FittingUtils *fitter = nullptr;
  FitResultDetailed result;

  Float_t fit_low, fit_high;

  if (peak_name == "Am_59keV") {
    fit_low = 2000;
    fit_high = 9000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kTRUE, kTRUE,
                              kTRUE, kTRUE);
  }

  result = fitter->FitPeakDetailed(input_name, peak_name);
  delete hist;
  delete fitter;
  return result;
}

FitResultStandard FitSinglePeak(const TString input_name,
                                const TString peak_name,
                                const Float_t expected_mu) {

  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << input_name << ".root" << std::endl;
    return {};
  }

  TH1F *hist = static_cast<TH1F *>(file->Get("long_integral"));
  if (!hist) {
    std::cerr << "Cannot find 'Pulse Integral' histogram in " << input_name
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

  Float_t fit_low, fit_high;

  if (peak_name == "Am_59keV") {
    fit_low = 2500;
    fit_high = 7000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Cs_662keV") {
    fit_low = 34000;
    fit_high = 45000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Na_511keV") {
    fit_low = 25000;
    fit_high = 36000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Na_1274keV") {
    fit_low = 72000;
    fit_high = 83000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
  } else if (peak_name == "Co_1332keV") {
    fit_low = 77000;
    fit_high = 83000;
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

  // Am-241 59.5 keV
  FitResultDetailed am_result =
      FitSinglePeakDetailed(Constants::AM241, "Am_59keV", 4000);
  cal_data.peak_names.push_back("Am_59keV");
  cal_data.mu.push_back(am_result.mu);
  cal_data.mu_errors.push_back(am_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_AM241_59KEV);
  cal_data.reduced_chi2.push_back(am_result.reduced_chi2);

  // Na-22 511 keV
  FitResultStandard na511_result =
      FitSinglePeak(Constants::NA22, "Na_511keV", 39000);
  cal_data.peak_names.push_back("Na_511keV");
  cal_data.mu.push_back(na511_result.mu);
  cal_data.mu_errors.push_back(na511_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_NA22_511KEV);
  cal_data.reduced_chi2.push_back(na511_result.reduced_chi2);

  // Cs-137 662 keV
  FitResultStandard cs_result =
      FitSinglePeak(Constants::CS137, "Cs_662keV", 30251);
  cal_data.peak_names.push_back("Cs_662keV");
  cal_data.mu.push_back(cs_result.mu);
  cal_data.mu_errors.push_back(cs_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_CS137_662KEV);
  cal_data.reduced_chi2.push_back(cs_result.reduced_chi2);

  // Na-22 1274.5 keV
  FitResultStandard na1274_result =
      FitSinglePeak(Constants::NA22, "Na_1274keV", 76000);
  cal_data.peak_names.push_back("Na_1274keV");
  cal_data.mu.push_back(na1274_result.mu);
  cal_data.mu_errors.push_back(na1274_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_NA22_1274KEV);
  cal_data.reduced_chi2.push_back(na1274_result.reduced_chi2);

  // Co-60 1332.5 keV
  FitResultStandard co_result =
      FitSinglePeak(Constants::CO60, "Co_1332keV", 78295);
  cal_data.peak_names.push_back("Co_1332keV");
  cal_data.mu.push_back(co_result.mu);
  cal_data.mu_errors.push_back(co_result.mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_CO60_1332KEV);
  cal_data.reduced_chi2.push_back(co_result.reduced_chi2);

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

  TF1 *calibration_fit = new TF1("linear", "pol1", 0, Constants::ADC_MAX);
  calibration_fit->SetParameter(0, 0);
  calibration_fit->SetParLimits(0, -1e-2, 1e-2);
  calibration_fit->SetParameter(1, 0.076);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");
  calibration_fit->Draw("SAME");
  PlottingUtils::SaveFigure(canvas, "calibration.png", kFALSE);

  delete canvas;
  return calibration_fit;
}

void LongIntegralToLightOutput(const std::vector<TString> &input_names,
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

    Float_t long_integral, light_output_keVee;
    features_tree->SetBranchAddress("long_integral", &long_integral);
    features_tree->Branch("light_output", &light_output_keVee,
                          "light_output/F");

    Int_t num_entries = features_tree->GetEntries();

    for (Int_t j = 0; j < num_entries; j++) {
      features_tree->GetEntry(j);
      light_output_keVee = calibration_function->Eval(long_integral);
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

  std::vector<TString> datasets_to_calibrate = Constants::ALL_OUTPUT_NAMES;
}
