#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSystem.h>
#include <cmath>
#include <iomanip>
#include <vector>

const Float_t E_AM241 = 59.5409;
const Float_t E_BA133_53 = 53.16;
const Float_t E_BA133_81 = 80.98;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;
const Float_t E_CD114M = 95.9023;

struct CalibrationData {
  std::vector<Float_t> mu;
  std::vector<Float_t> mu_errors;
  std::vector<Float_t> calibration_values_keV;
  std::vector<Float_t> reduced_chi2;
  std::vector<TString> run_names;
  std::vector<TString> dataset_names;
};

struct DriftCorrectionResult {
  Float_t correction_factor;
  Float_t correction_factor_error;
  std::vector<Float_t> individual_factors;
  std::vector<Float_t> individual_errors;
  std::vector<TString> peak_names;
};

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label);

FitResultDetailed FitSinglePeak(const TString input_name,
                                const TString peak_name,
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
  FitResultDetailed result;

  if (peak_name == "Am_59.5keV") {
    fitter =
        new FittingUtils(zoomedHist, 50, 65, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE);
  }
  if (peak_name == "Ba_80.98keV") {
    if (input_name == Constants::POSTREACTOR_BA133_20260115) {
      fitter = new FittingUtils(zoomedHist, 72, 87, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else {
      fitter = new FittingUtils(zoomedHist, 70, 90, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    }
  }
  if (peak_name == "Cd114m_95.9keV") {
    if (input_name == Constants::NOSHIELDBACKGROUND_5PERCENT_20260115) {
      fitter = new FittingUtils(zoomedHist, 91, 100, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else {
      fitter = new FittingUtils(zoomedHist, 91, 103, kFALSE, kTRUE, kTRUE,
                                kTRUE, kTRUE);
    }
  }
  result = fitter->FitPeakDetailed(input_name, peak_name);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResultDoublePeakDetailed
FitDoublePeak(const TString input_name, const TString peak_name,
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
  FitResultDoublePeakDetailed result;

  if (peak_name == "Am_59.5_Ba_53.16keV") {
    fitter = new FittingUtils(zoomedHist, 45, 70, kFALSE, kTRUE, kTRUE, kTRUE,
                              kTRUE);
  }
  if (peak_name == "Pb_KAlpha") {
    if (input_name == Constants::CDSHIELDBACKGROUND_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 66, 81, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260114) {
      fitter = new FittingUtils(zoomedHist, 66, 82, kTRUE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else if (input_name == Constants::CUSHIELDBACKGROUND_10PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 65, 82, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result =
      fitter->FitDoublePeakDetailed(input_name, peak_name, mu1_init, mu2_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

FitResultTriplePeakDetailed
FitTriplePeak(const TString input_name, const TString peak_name,
              const FitResultDoublePeakDetailed &constrained_peaks,
              const Float_t mu3_init, const Bool_t use_calibrated = kFALSE) {

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
  FitResultTriplePeakDetailed result;

  if (peak_name == "Ge") {
    if (input_name == Constants::CDSHIELDSIGNAL_25PERCENT_20260113) {
      fitter = new FittingUtils(zoomedHist, 63, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
    } else if (input_name == Constants::CUSHIELDSIGNAL_10PERCENT_20260114) {
      fitter = new FittingUtils(zoomedHist, 60, 80, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);

    } else
      fitter = new FittingUtils(zoomedHist, 65, 81, kFALSE, kTRUE, kTRUE, kTRUE,
                                kTRUE);
  }

  result = fitter->FitTriplePeakDetailed(input_name, peak_name,
                                         constrained_peaks, mu3_init);
  delete zoomedHist;
  delete fitter;
  return result;
}

// Calculate multiplicative drift correction factor by comparing measured peak
// positions to reference positions. Returns weighted average of
// (reference/measured) ratios.
DriftCorrectionResult
CalculateDriftCorrection(const std::vector<Float_t> &ref_positions,
                         const std::vector<Float_t> &ref_errors,
                         const std::vector<Float_t> &measured_positions,
                         const std::vector<Float_t> &measured_errors,
                         const std::vector<TString> &peak_names) {

  DriftCorrectionResult result;
  result.peak_names = peak_names;

  Float_t weighted_sum = 0;
  Float_t weight_sum = 0;

  for (size_t i = 0; i < ref_positions.size(); ++i) {
    if (measured_positions[i] > 0 && ref_positions[i] > 0) {
      Float_t factor = ref_positions[i] / measured_positions[i];
      // Error propagation for ratio
      Float_t rel_err_ref = ref_errors[i] / ref_positions[i];
      Float_t rel_err_meas = measured_errors[i] / measured_positions[i];
      Float_t factor_error = factor * std::sqrt(rel_err_ref * rel_err_ref +
                                                rel_err_meas * rel_err_meas);

      result.individual_factors.push_back(factor);
      result.individual_errors.push_back(factor_error);

      Float_t weight = 1.0 / (factor_error * factor_error);
      weighted_sum += factor * weight;
      weight_sum += weight;

      std::cout << "  " << peak_names[i] << ": factor = " << std::fixed
                << std::setprecision(6) << factor << " +/- " << factor_error
                << " (mu = " << measured_positions[i] << " +/- "
                << measured_errors[i] << " keV)" << std::endl;
    }
  }

  if (weight_sum > 0) {
    result.correction_factor = weighted_sum / weight_sum;
    result.correction_factor_error = std::sqrt(1.0 / weight_sum);
  } else {
    result.correction_factor = 1.0;
    result.correction_factor_error = 0.0;
  }

  std::cout << "  Weighted average correction factor: " << std::fixed
            << std::setprecision(6) << result.correction_factor << " +/- "
            << result.correction_factor_error << std::endl;

  return result;
}

TF1 *CreateAndSaveCalibration(std::vector<Float_t> mu,
                              std::vector<Float_t> calibration_values_keV,
                              std::vector<Float_t> mu_errors,
                              TString date_label) {

  Int_t size = calibration_values_keV.size();

  TGraph *calibration_curve =
      new TGraph(size, mu.data(), calibration_values_keV.data());
  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(
      calibration_curve, kBlue,
      "; Precalibrated Energy [keV]; Deposited Energy [keV]");

  calibration_curve->GetXaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetYaxis()->SetRangeUser(-5, 100);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->SetMarkerStyle(5);
  calibration_curve->SetMarkerSize(2);
  calibration_curve->Draw("AP");

  TF1 *calibration_fit = new TF1("linear_" + date_label, "pol1", -10, 100);
  calibration_fit->FixParameter(0, 0);
  calibration_fit->SetParameter(1, 1);
  calibration_fit->SetNpx(1000);

  TFitResultPtr fit_result = calibration_curve->Fit(calibration_fit, "LRE");

  calibration_fit->Draw("SAME");

  PlottingUtils::SaveFigure(canvas, "calibration_" + date_label + ".png",
                            kFALSE);
  return calibration_fit;
}

void PulseHeightToDepositedEnergy(std::vector<TString> input_names,
                                  TF1 *calibration_function,
                                  TString date_label) {
  Int_t entries = input_names.size();

  TString calibration_function_filepath =
      "root_files/calibration_function_" + date_label + ".root";
  TFile *calibration_file =
      new TFile(calibration_function_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);

  for (Int_t i = 0; i < entries; i++) {
    TString input_name = input_names[i];

    TH1F *hist = new TH1F(PlottingUtils::GetRandomName(),
                          Form("; Deposited Energy [keV]; Counts / %d eV",
                               Constants::BIN_WIDTH_EV),
                          Constants::HIST_NBINS, Constants::HIST_XMIN,
                          Constants::HIST_XMAX);

    TH1F *zoomedHist = new TH1F(PlottingUtils::GetRandomName(),
                                Form("; Deposited Energy [keV]; Counts / %d eV",
                                     Constants::BIN_WIDTH_EV),
                                Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                                Constants::ZOOMED_XMAX);

    TH1F *peakHist = new TH1F(PlottingUtils::GetRandomName(),
                              Form("; Deposited Energy [keV]; Counts / %d eV",
                                   Constants::BIN_WIDTH_EV),
                              Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                              Constants::PEAK_XMAX);

    TString output_filepath = "root_files/" + input_name + ".root";

    TFile *output = new TFile(output_filepath, "UPDATE");

    output->cd();

    Bool_t isFiltered = input_name.Contains("_filtered");

    TTree *tree = nullptr;
    TString treeName;
    TString energyBranchName;

    if (isFiltered) {
      treeName = "bef_tree";
      energyBranchName = "energykeV";
      tree = static_cast<TTree *>(output->Get(treeName));
    } else {
      treeName = "bef_tree_event_summary";
      energyBranchName = "totalEnergykeV";
      tree = static_cast<TTree *>(output->Get(treeName));
    }

    if (!tree) {
      std::cerr << "ERROR: Could not find tree '" << treeName << "' in file "
                << output_filepath << std::endl;
      output->Close();
      delete hist;
      delete zoomedHist;
      delete peakHist;
      continue;
    }

    Float_t energy, deposited_energy_keV;
    tree->SetBranchAddress(energyBranchName, &energy);
    tree->Branch("deposited_energy", &deposited_energy_keV,
                 "deposited_energy/F");

    Int_t num_entries = tree->GetEntries();

    tree->LoadBaskets();
    for (Int_t j = 0; j < num_entries; j++) {
      tree->GetEntry(j);
      deposited_energy_keV = calibration_function->Eval(energy);
      tree->GetBranch("deposited_energy")->Fill();
      hist->Fill(deposited_energy_keV);
      zoomedHist->Fill(deposited_energy_keV);
      peakHist->Fill(deposited_energy_keV);
    }

    PlottingUtils::ConfigureHistogram(hist, kViolet);

    tree->Write("", TObject::kOverwrite);
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);
    output->Close();
    delete hist;
    delete zoomedHist;
    delete peakHist;
  }
  calibration_file->Close();
}

struct ReferenceCalibration {
  TF1 *calibration_function;
  Float_t ref_pb_ka1;
  Float_t ref_pb_ka1_error;
  Float_t ref_pb_ka2;
  Float_t ref_pb_ka2_error;
  Float_t ref_cd114m;
  Float_t ref_cd114m_error;
};

ReferenceCalibration GetReferenceCalibration_20260113() {
  std::cout << "Building Reference Calibration for Jan 13" << std::endl;

  CalibrationData cal_data;

  cal_data.run_names.push_back("Zero Point");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);

  FitResultDetailed am241_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_20260113, "Am_59.5keV", E_AM241);
  cal_data.run_names.push_back("Post-reactor Am-241");
  cal_data.mu.push_back(am241_result.mu);
  cal_data.mu_errors.push_back(am241_result.mu_error);
  cal_data.calibration_values_keV.push_back(E_AM241);
  cal_data.reduced_chi2.push_back(am241_result.reduced_chi2);

  FitResultDoublePeakDetailed cu_bkg_pb =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                    "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
  cal_data.run_names.push_back("Cu Shield Bkg Pb-Ka1 (reference)");
  cal_data.mu.push_back(cu_bkg_pb.peak1.mu);
  cal_data.mu_errors.push_back(cu_bkg_pb.peak1.mu_error);
  cal_data.calibration_values_keV.push_back(E_PB_KA1);
  cal_data.reduced_chi2.push_back(cu_bkg_pb.reduced_chi2);

  cal_data.run_names.push_back("Cu Shield Bkg Pb-Ka2 (reference)");
  cal_data.mu.push_back(cu_bkg_pb.peak2.mu);
  cal_data.mu_errors.push_back(cu_bkg_pb.peak2.mu_error);
  cal_data.calibration_values_keV.push_back(E_PB_KA2);
  cal_data.reduced_chi2.push_back(-1);

  FitResultDetailed cu_bkg_cd =
      FitSinglePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                    "Cd114m_95.9keV", E_CD114M);
  cal_data.run_names.push_back("Cu Shield Bkg Cd-114m (reference)");
  cal_data.mu.push_back(cu_bkg_cd.mu);
  cal_data.mu_errors.push_back(cu_bkg_cd.mu_error);
  cal_data.calibration_values_keV.push_back(E_CD114M);
  cal_data.reduced_chi2.push_back(cu_bkg_cd.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260113_reference");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260113_reference");

  ReferenceCalibration ref;
  ref.calibration_function = cal_func;
  ref.ref_pb_ka1 = cu_bkg_pb.peak1.mu;
  ref.ref_pb_ka1_error = cu_bkg_pb.peak1.mu_error;
  ref.ref_pb_ka2 = cu_bkg_pb.peak2.mu;
  ref.ref_pb_ka2_error = cu_bkg_pb.peak2.mu_error;
  ref.ref_cd114m = cu_bkg_cd.mu;
  ref.ref_cd114m_error = cu_bkg_cd.mu_error;

  return ref;
}

void ApplyCalibration_CuShield_20260113(TF1 *calibration_function) {
  std::cout << "Applying calibration to CuShield runs + PostReactor "
               "(Jan 13)"
            << std::endl;

  std::vector<TString> datasets = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
      Constants::POSTREACTOR_AM241_20260113};

  PulseHeightToDepositedEnergy(datasets, calibration_function,
                               "20260113_cushield");
}

void ApplyCalibration_CdShield_20260113(const ReferenceCalibration &ref) {
  std::cout << "Applying drift-corrected calibration to CdShield runs "
               "(Jan 13)"
            << std::endl;

  std::vector<std::pair<TString, TString>> cd_configs = {
      {Constants::CDSHIELDBACKGROUND_10PERCENT_20260113, "CdShield_10%"},
      {Constants::CDSHIELDBACKGROUND_25PERCENT_20260113, "CdShield_25%"}};

  for (const auto &config : cd_configs) {
    TString bkg_name = config.first;
    TString label = config.second;

    std::cout << "Processing " << label << ":" << std::endl;

    FitResultDoublePeakDetailed pb_fit =
        FitDoublePeak(bkg_name, "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
    FitResultDetailed cd_fit =
        FitSinglePeak(bkg_name, "Cd114m_95.9keV", E_CD114M);

    std::vector<Float_t> ref_pos = {ref.ref_pb_ka1, ref.ref_pb_ka2,
                                    ref.ref_cd114m};
    std::vector<Float_t> ref_err = {ref.ref_pb_ka1_error, ref.ref_pb_ka2_error,
                                    ref.ref_cd114m_error};
    std::vector<Float_t> meas_pos = {pb_fit.peak1.mu, pb_fit.peak2.mu,
                                     cd_fit.mu};
    std::vector<Float_t> meas_err = {pb_fit.peak1.mu_error,
                                     pb_fit.peak2.mu_error, cd_fit.mu_error};
    std::vector<TString> names = {"Pb-Ka1", "Pb-Ka2", "Cd-114m"};

    DriftCorrectionResult drift =
        CalculateDriftCorrection(ref_pos, ref_err, meas_pos, meas_err, names);

    Float_t original_slope = ref.calibration_function->GetParameter(1);
    Float_t corrected_slope = original_slope * drift.correction_factor;

    TString func_name = "linear_20260113_" + label;
    TF1 *corrected_cal = new TF1(func_name, "pol1", -10, 100);
    corrected_cal->FixParameter(0, 0);
    corrected_cal->SetParameter(1, corrected_slope);

    std::cout << "  Original slope: " << std::fixed << std::setprecision(6)
              << original_slope << ", Corrected slope: " << corrected_slope
              << std::endl;

    TString signal_name;
    if (bkg_name.Contains("10Percent")) {
      signal_name = Constants::CDSHIELDSIGNAL_10PERCENT_20260113;
    } else {
      signal_name = Constants::CDSHIELDSIGNAL_25PERCENT_20260113;
    }

    std::vector<TString> datasets = {bkg_name, signal_name};
    TString date_label = "20260113_" + label;
    date_label.ReplaceAll("%", "pct");
    PulseHeightToDepositedEnergy(datasets, corrected_cal, date_label);
  }
}

void ApplyCalibration_20260114(const ReferenceCalibration &ref_jan13) {
  std::cout << "Applying drift-corrected calibration to Jan 14 runs"
            << std::endl;

  FitResultDoublePeakDetailed pb_fit =
      FitDoublePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                    "Pb_KAlpha", E_PB_KA1, E_PB_KA2);
  FitResultDetailed cd_fit =
      FitSinglePeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                    "Cd114m_95.9keV", E_CD114M);

  std::vector<Float_t> ref_pos = {ref_jan13.ref_pb_ka1, ref_jan13.ref_pb_ka2,
                                  ref_jan13.ref_cd114m};
  std::vector<Float_t> ref_err = {ref_jan13.ref_pb_ka1_error,
                                  ref_jan13.ref_pb_ka2_error,
                                  ref_jan13.ref_cd114m_error};
  std::vector<Float_t> meas_pos = {pb_fit.peak1.mu, pb_fit.peak2.mu, cd_fit.mu};
  std::vector<Float_t> meas_err = {pb_fit.peak1.mu_error, pb_fit.peak2.mu_error,
                                   cd_fit.mu_error};
  std::vector<TString> names = {"Pb-Ka1", "Pb-Ka2", "Cd-114m"};

  DriftCorrectionResult drift =
      CalculateDriftCorrection(ref_pos, ref_err, meas_pos, meas_err, names);

  Float_t original_slope = ref_jan13.calibration_function->GetParameter(1);
  Float_t corrected_slope = original_slope * drift.correction_factor;

  TF1 *corrected_cal = new TF1("linear_20260114", "pol1", -10, 100);
  corrected_cal->FixParameter(0, 0);
  corrected_cal->SetParameter(1, corrected_slope);

  std::cout << "  Original slope (Jan 13): " << std::fixed
            << std::setprecision(6) << original_slope
            << ", Corrected slope (Jan 14): " << corrected_slope << std::endl;

  std::vector<TString> datasets = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_90PERCENT_20260114};

  PulseHeightToDepositedEnergy(datasets, corrected_cal, "20260114");
}

ReferenceCalibration GetReferenceCalibration_20260115() {
  std::cout << "Building Reference Calibration for Jan 15" << std::endl;

  CalibrationData cal_data;

  cal_data.run_names.push_back("Zero Point");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);

  FitResultDetailed am241_result = FitSinglePeak(
      Constants::POSTREACTOR_AM241_20260115, "Am_59.5keV", E_AM241);
  cal_data.run_names.push_back("Post-reactor Am-241");
  cal_data.mu.push_back(am241_result.mu);
  cal_data.mu_errors.push_back(am241_result.mu_error);
  cal_data.calibration_values_keV.push_back(E_AM241);
  cal_data.reduced_chi2.push_back(am241_result.reduced_chi2);

  FitResultDetailed ba133_result = FitSinglePeak(
      Constants::POSTREACTOR_BA133_20260115, "Ba_80.98keV", E_BA133_81);
  cal_data.run_names.push_back("Post-reactor Ba-133 80.98keV");
  cal_data.mu.push_back(ba133_result.mu);
  cal_data.mu_errors.push_back(ba133_result.mu_error);
  cal_data.calibration_values_keV.push_back(E_BA133_81);
  cal_data.reduced_chi2.push_back(ba133_result.reduced_chi2);

  FitResultDetailed cd_fit =
      FitSinglePeak(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                    "Cd114m_95.9keV", E_CD114M);
  cal_data.run_names.push_back("No Shield Bkg Cd-114m (reference)");
  cal_data.mu.push_back(cd_fit.mu);
  cal_data.mu_errors.push_back(cd_fit.mu_error);
  cal_data.calibration_values_keV.push_back(E_CD114M);
  cal_data.reduced_chi2.push_back(cd_fit.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260115_reference");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260115_reference");

  ReferenceCalibration ref;
  ref.calibration_function = cal_func;
  ref.ref_pb_ka1 = 0; // No Pb on Jan 15
  ref.ref_pb_ka1_error = 0;
  ref.ref_pb_ka2 = 0;
  ref.ref_pb_ka2_error = 0;
  ref.ref_cd114m = cd_fit.mu;
  ref.ref_cd114m_error = cd_fit.mu_error;

  return ref;
}

void ApplyCalibration_20260115(TF1 *calibration_function) {
  std::cout << "Applying calibration to Jan 15 runs" << std::endl;

  std::vector<TString> datasets = {
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
      Constants::POSTREACTOR_AM241_20260115,
      Constants::POSTREACTOR_BA133_20260115, Constants::SHUTTERCLOSED_20260115};

  PulseHeightToDepositedEnergy(datasets, calibration_function, "20260115");
}

void ApplyCalibration_20260116(const ReferenceCalibration &ref_jan15) {
  std::cout << "Applying drift-corrected calibration to Jan 16 runs"
            << std::endl;

  FitResultDetailed cd_fit = FitSinglePeak(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      "Cd114m_95.9keV", E_CD114M);

  std::vector<Float_t> ref_pos = {ref_jan15.ref_cd114m};
  std::vector<Float_t> ref_err = {ref_jan15.ref_cd114m_error};
  std::vector<Float_t> meas_pos = {cd_fit.mu};
  std::vector<Float_t> meas_err = {cd_fit.mu_error};
  std::vector<TString> names = {"Cd-114m"};

  DriftCorrectionResult drift =
      CalculateDriftCorrection(ref_pos, ref_err, meas_pos, meas_err, names);

  Float_t original_slope = ref_jan15.calibration_function->GetParameter(1);
  Float_t corrected_slope = original_slope * drift.correction_factor;

  TF1 *corrected_cal = new TF1("linear_20260116", "pol1", -10, 100);
  corrected_cal->FixParameter(0, 0);
  corrected_cal->SetParameter(1, corrected_slope);

  std::cout << "  Original slope (Jan 15): " << std::fixed
            << std::setprecision(6) << original_slope
            << ", Corrected slope (Jan 16): " << corrected_slope << std::endl;

  std::vector<TString> datasets = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};

  PulseHeightToDepositedEnergy(datasets, corrected_cal, "20260116");
}

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label) {
  std::cout << "Calibration Points for " << date_label << std::endl;
  for (size_t i = 0; i < cal_data.mu.size(); ++i) {
    std::cout << std::left << std::setw(45) << cal_data.run_names[i] << ": "
              << std::fixed << std::setprecision(4) << cal_data.mu[i] << " +/- "
              << cal_data.mu_errors[i] << " keV";
    if (cal_data.reduced_chi2[i] > 0) {
      std::cout << " (chi2/ndf = " << std::setprecision(3)
                << cal_data.reduced_chi2[i] << ")";
      std::cout << std::endl;
    }
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences();
  std::cout << "Jan 13 Calibration" << std::endl;
  ReferenceCalibration ref_jan13 = GetReferenceCalibration_20260113();
  ApplyCalibration_CuShield_20260113(ref_jan13.calibration_function);
  ApplyCalibration_CdShield_20260113(ref_jan13);

  std::cout << "Jan 14 Calibration" << std::endl;
  ApplyCalibration_20260114(ref_jan13);

  std::cout << "Jan 15 Calibration" << std::endl;
  ReferenceCalibration ref_jan15 = GetReferenceCalibration_20260115();
  ApplyCalibration_20260115(ref_jan15.calibration_function);

  std::cout << "Jan 16 Calibration" << std::endl;
  ApplyCalibration_20260116(ref_jan15);
}
