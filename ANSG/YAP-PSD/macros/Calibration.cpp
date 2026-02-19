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

FitResult FitSinglePeak(const TString input_name, const TString peak_name,
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
  FitResult result;

  Float_t fit_low, fit_high;
  if (peak_name == "Am_59keV") {
    fit_low = 500;
    fit_high = 9500;
    fitter = new FittingUtils(hist, fit_low, fit_high);
  } else if (peak_name == "Na_511keV") {
    fit_low = 26000;
    fit_high = 35000;
    fitter = new FittingUtils(hist, fit_low, fit_high);
  } else if (peak_name == "Cs_662keV") {
    fit_low = 35500;
    fit_high = 47000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kTRUE, kTRUE);
  } else if (peak_name == "Na_1274keV") {
    fit_low = 70000;
    fit_high = 84000;
    fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kTRUE);
  }

  result = fitter->FitPeak(input_name, peak_name);
  delete hist;
  delete fitter;
  return result;
}

CalibrationData FitCalibrationPeaks() {
  CalibrationData cal_data;

  cal_data.peak_names.push_back("Zero");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);

  FitResult am_result = FitSinglePeak(Constants::AM241, "Am_59keV", 3700);
  cal_data.peak_names.push_back("Am_59keV");
  cal_data.mu.push_back(am_result.peaks.at(0).mu);
  cal_data.mu_errors.push_back(am_result.peaks.at(0).mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_AM241_59KEV);
  cal_data.reduced_chi2.push_back(am_result.reduced_chi2);

  FitResult na511_result = FitSinglePeak(Constants::NA22, "Na_511keV", 31000);
  cal_data.peak_names.push_back("Na_511keV");
  cal_data.mu.push_back(na511_result.peaks.at(0).mu);
  cal_data.mu_errors.push_back(na511_result.peaks.at(0).mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_NA22_511KEV);
  cal_data.reduced_chi2.push_back(na511_result.reduced_chi2);

  FitResult cs_result = FitSinglePeak(Constants::CS137, "Cs_662keV", 36850);
  cal_data.peak_names.push_back("Cs_662keV");
  cal_data.mu.push_back(cs_result.peaks.at(0).mu);
  cal_data.mu_errors.push_back(cs_result.peaks.at(0).mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_CS137_662KEV);
  cal_data.reduced_chi2.push_back(cs_result.reduced_chi2);

  FitResult na1274_result = FitSinglePeak(Constants::NA22, "Na_1274keV", 70000);
  cal_data.peak_names.push_back("Na_1274keV");
  cal_data.mu.push_back(na1274_result.peaks.at(0).mu);
  cal_data.mu_errors.push_back(na1274_result.peaks.at(0).mu_error);
  cal_data.calibration_values_keV.push_back(Constants::E_NA22_1274KEV);
  cal_data.reduced_chi2.push_back(na1274_result.reduced_chi2);

  return cal_data;
}

void PrintCalibrationSummary(const CalibrationData &cal_data) {
  for (size_t i = 0; i < cal_data.mu.size(); ++i) {
    std::cout << cal_data.peak_names[i] << ": " << std::fixed
              << std::setprecision(2) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " a.u. -> "
              << cal_data.calibration_values_keV[i] << " keV" << std::endl;
  }
}

TF1 *CreateAndSaveCalibration(const CalibrationData &cal_data) {

  Int_t size = cal_data.calibration_values_keV.size();

  TGraphErrors *calibration_curve = new TGraphErrors(
      size, cal_data.mu.data(), cal_data.calibration_values_keV.data(),
      cal_data.mu_errors.data(), nullptr);

  TF1 *quadratic_fit = new TF1("calibration_function", "pol2", -10, 90000);
  quadratic_fit->SetParameter(0, 0);
  quadratic_fit->SetParameter(1, 0.02);

  calibration_curve->Fit(quadratic_fit, "LRE");
  quadratic_fit->SetNpx(10000);

  std::vector<Float_t> mu_no_am, mu_err_no_am, energy_no_am;
  for (Int_t i = 0; i < size; ++i) {
    if (cal_data.peak_names[i] == "Am_59keV")
      continue;
    mu_no_am.push_back(cal_data.mu[i]);
    mu_err_no_am.push_back(cal_data.mu_errors[i]);
    energy_no_am.push_back(cal_data.calibration_values_keV[i]);
  }

  TGraphErrors *graph_no_am =
      new TGraphErrors(mu_no_am.size(), mu_no_am.data(), energy_no_am.data(),
                       mu_err_no_am.data(), nullptr);

  TF1 *linear_fit = new TF1("linear_fit", "pol1", -10, 90000);
  linear_fit->SetParameter(0, 0);
  linear_fit->SetParameter(1, 0.02);
  graph_no_am->Fit(linear_fit, "LRE");
  linear_fit->SetNpx(10000);

  TCanvas *canvas = new TCanvas("c_cal_curve", "Calibration Curve", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);

  calibration_curve->SetMarkerStyle(20);
  calibration_curve->SetMarkerSize(1.2);
  calibration_curve->SetMarkerColor(kBlack);
  calibration_curve->SetLineColor(kBlack);
  calibration_curve->SetTitle("");
  calibration_curve->GetXaxis()->SetTitle("Pulse Integral [a.u.]");
  calibration_curve->GetYaxis()->SetTitle("Deposited Energy [keV]");
  calibration_curve->GetXaxis()->SetTitleSize(0.06);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->GetYaxis()->SetTitleSize(0.06);
  calibration_curve->GetXaxis()->SetLabelSize(0.06);
  calibration_curve->GetYaxis()->SetLabelSize(0.06);
  calibration_curve->GetXaxis()->SetTitleOffset(1.2);
  calibration_curve->GetYaxis()->SetTitleOffset(1.2);

  Double_t x_min, x_max, y_min, y_max;
  calibration_curve->ComputeRange(x_min, y_min, x_max, y_max);
  calibration_curve->GetYaxis()->SetRangeUser(-10, y_max * 1.1);
  calibration_curve->GetXaxis()->SetRangeUser(-10, x_max * 1.1);

  quadratic_fit->GetXaxis()->SetRangeUser(-10, x_max * 1.1);
  quadratic_fit->SetLineColor(kRed);
  quadratic_fit->SetLineStyle(1);

  linear_fit->GetXaxis()->SetRangeUser(-10, x_max * 1.1);
  linear_fit->SetLineColor(kBlue);
  linear_fit->SetLineStyle(1);

  calibration_curve->GetListOfFunctions()->Clear();

  calibration_curve->Draw("APE"); // points only now, establishes axes
  quadratic_fit->Draw("SAME");
  linear_fit->Draw("SAME");
  calibration_curve->Draw("PE SAME"); // points on top

  TLegend *leg = new TLegend(0.2, 0.65, 0.6, 0.84);
  leg->SetMargin(0.2);
  leg->AddEntry(calibration_curve, "Calibration points", "pe");
  leg->AddEntry(quadratic_fit, "Quadratic fit", "l");
  leg->AddEntry(linear_fit, "Linear fit (high energy)", "l");
  leg->Draw();

  PlottingUtils::SaveFigure(canvas, "calibration.png",
                            PlotSaveOptions::kLINEAR);

  delete graph_no_am;
  delete canvas;
  return quadratic_fit;
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
    PlottingUtils::SaveFigure(canvas, input_name + "_light_output.png",
                              PlotSaveOptions::kLOG);

    features_tree->Write("", TObject::kOverwrite);
    light_output_hist->Write("Light Output", TObject::kOverwrite);
    output->Close();
    delete output;
    delete canvas;
    delete light_output_hist;
  }
}

void Calibration() {
  Bool_t recalibrate = kTRUE;
  InitUtils::SetROOTPreferences();

  CalibrationData cal_data = FitCalibrationPeaks();

  PrintCalibrationSummary(cal_data);

  TF1 *calibration_function = CreateAndSaveCalibration(cal_data);

  if (recalibrate)
    LongIntegralToLightOutput(Constants::ALL_OUTPUT_NAMES,
                              calibration_function);
}
