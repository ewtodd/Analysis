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
};

void PrintCalibrationSummary(const CalibrationData &cal_data,
                             TString date_label);

TH1F *LoadHistogram(const TString input_name,
                    const Bool_t use_calibrated = kFALSE) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }
  std::cout << "FOUND " << input_name << ".root" << std::endl;

  TString histName = use_calibrated ? "calibrated_zoomedHist" : "zoomedHist";
  TH1F *hist = static_cast<TH1F *>(file->Get(histName));
  hist->SetDirectory(0);
  file->Close();
  delete file;
  return hist;
}

FitResult FitCalibrationPeak(const TString input_name, const TString peak_name,
                             const Bool_t use_calibrated,
                             const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, use_calibrated);
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

FitResult FitPbKAlpha(const TString input_name, const Bool_t use_calibrated,
                      const Bool_t interactive) {
  TH1F *hist = LoadHistogram(input_name, use_calibrated);
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

TF1 *CreateAndSaveCalibration(
    const std::vector<Float_t> &mu,
    const std::vector<Float_t> &calibration_values_keV,
    const std::vector<Float_t> &mu_errors, const TString date_label) {

  Int_t size = calibration_values_keV.size();

  TGraph *calibration_curve =
      new TGraph(size, mu.data(), calibration_values_keV.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas();
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

  PlottingUtils::SaveFigure(canvas, "calibration_" + date_label, "",
                            PlotSaveOptions::kLINEAR);
  return calibration_fit;
}

void AddZeroPoint(CalibrationData &cal_data) {
  cal_data.run_names.push_back("Zero Point");
  cal_data.mu.push_back(0);
  cal_data.mu_errors.push_back(0);
  cal_data.calibration_values_keV.push_back(0);
  cal_data.reduced_chi2.push_back(0);
}

void AddCalibrationPoint(CalibrationData &cal_data, const TString run_name,
                         const Float_t mu, const Float_t mu_error,
                         const Float_t true_energy,
                         const Float_t reduced_chi2) {
  cal_data.run_names.push_back(run_name);
  cal_data.mu.push_back(mu);
  cal_data.mu_errors.push_back(mu_error);
  cal_data.calibration_values_keV.push_back(true_energy);
  cal_data.reduced_chi2.push_back(reduced_chi2);
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
    }
    std::cout << std::endl;
  }
}

TF1 *BuildCalibration_20260113(const Bool_t use_calibrated,
                               const Bool_t interactive) {
  std::cout << "Building Calibration for Jan 13" << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_20260113, "Am_59.5keV",
                         use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult cu_bkg_pb =
      FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                  use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                      cu_bkg_pb.peaks.at(0).mu, cu_bkg_pb.peaks.at(0).mu_error,
                      E_PB_KA1, cu_bkg_pb.reduced_chi2);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                      cu_bkg_pb.peaks.at(1).mu, cu_bkg_pb.peaks.at(1).mu_error,
                      E_PB_KA2, -1);

  FitResult cu_bkg_cd =
      FitCalibrationPeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
                         "Cd114m_95.9keV", use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Cd-114m",
                      cu_bkg_cd.peaks.at(0).mu, cu_bkg_cd.peaks.at(0).mu_error,
                      E_CD114M, cu_bkg_cd.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260113");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260113");
  return cal_func;
}

TF1 *BuildCalibration_20260114(const Bool_t use_calibrated,
                               const Bool_t interactive) {
  std::cout << "Building Calibration for Jan 14" << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult cu_bkg_pb =
      FitPbKAlpha(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                  use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka1",
                      cu_bkg_pb.peaks.at(0).mu, cu_bkg_pb.peaks.at(0).mu_error,
                      E_PB_KA1, cu_bkg_pb.reduced_chi2);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Pb-Ka2",
                      cu_bkg_pb.peaks.at(1).mu, cu_bkg_pb.peaks.at(1).mu_error,
                      E_PB_KA2, -1);

  FitResult cu_bkg_cd =
      FitCalibrationPeak(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
                         "Cd114m_95.9keV", use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "Cu Shield Bkg Cd-114m",
                      cu_bkg_cd.peaks.at(0).mu, cu_bkg_cd.peaks.at(0).mu_error,
                      E_CD114M, cu_bkg_cd.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260114");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260114");
  return cal_func;
}

TF1 *BuildCalibration_20260115(const Bool_t use_calibrated,
                               const Bool_t interactive) {
  std::cout << "Building Calibration for Jan 15" << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result =
      FitCalibrationPeak(Constants::POSTREACTOR_AM241_20260115, "Am_59.5keV",
                         use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba133_result =
      FitCalibrationPeak(Constants::POSTREACTOR_BA133_20260115, "Ba_80.98keV",
                         use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Ba-133 80.98keV", ba133_result.peaks.at(0).mu,
      ba133_result.peaks.at(0).mu_error, E_BA133_81, ba133_result.reduced_chi2);

  FitResult cd_fit =
      FitCalibrationPeak(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
                         "Cd114m_95.9keV", use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "No Shield Bkg Cd-114m", cd_fit.peaks.at(0).mu,
                      cd_fit.peaks.at(0).mu_error, E_CD114M,
                      cd_fit.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260115");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260115");
  return cal_func;
}

TF1 *BuildCalibration_20260116(const Bool_t use_calibrated,
                               const Bool_t interactive) {
  std::cout << "Building Calibration for Jan 16" << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult am241_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_BA133_20260116, "Am_59.5keV",
      use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Am-241", am241_result.peaks.at(0).mu,
      am241_result.peaks.at(0).mu_error, E_AM241, am241_result.reduced_chi2);

  FitResult ba133_result = FitCalibrationPeak(
      Constants::POSTREACTOR_AM241_BA133_20260116, "Ba_80.98keV",
      use_calibrated, interactive);
  AddCalibrationPoint(
      cal_data, "Post-reactor Ba-133 80.98keV", ba133_result.peaks.at(0).mu,
      ba133_result.peaks.at(0).mu_error, E_BA133_81,
      ba133_result.reduced_chi2);

  FitResult cd_fit = FitCalibrationPeak(
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      "Cd114m_95.9keV", use_calibrated, interactive);
  AddCalibrationPoint(cal_data, "GraphiteCastle Bkg Cd-114m",
                      cd_fit.peaks.at(0).mu, cd_fit.peaks.at(0).mu_error,
                      E_CD114M, cd_fit.reduced_chi2);

  PrintCalibrationSummary(cal_data, "20260116");

  TF1 *cal_func =
      CreateAndSaveCalibration(cal_data.mu, cal_data.calibration_values_keV,
                               cal_data.mu_errors, "20260116");
  return cal_func;
}

void PulseHeightToDepositedEnergy(const std::vector<TString> &input_names,
                                  TF1 *calibration_function,
                                  const TString date_label) {
  TString cal_filepath =
      "root_files/calibration_function_" + date_label + ".root";
  TFile *cal_file = new TFile(cal_filepath, "RECREATE");
  calibration_function->Write("calibration", TObject::kOverwrite);
  cal_file->Close();
  delete cal_file;

  Int_t n_files = input_names.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_name = input_names[i];
    TString filepath = "root_files/" + input_name + ".root";
    TFile *file = new TFile(filepath, "UPDATE");

    TTree *tree = static_cast<TTree *>(file->Get("bef_tree"));
    if (!tree) {
      std::cerr << "ERROR: Could not find bef_tree in " << filepath
                << std::endl;
      file->Close();
      delete file;
      continue;
    }

    Float_t energy;
    Float_t x, y, z;
    Int_t nInteractions;
    tree->SetBranchAddress("energykeV", &energy);
    tree->SetBranchAddress("xum", &x);
    tree->SetBranchAddress("yum", &y);
    tree->SetBranchAddress("zum", &z);
    tree->SetBranchAddress("nInteractions", &nInteractions);

    Float_t deposited_energy;
    TTree *cal_tree =
        new TTree("calibrated_tree", "Filtered and calibrated events");
    cal_tree->SetDirectory(file);
    cal_tree->Branch("deposited_energy", &deposited_energy,
                     "deposited_energy/F");

    TH1F *hist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Deposited Energy [keV]; Counts / %d eV",
             Constants::BIN_WIDTH_EV),
        Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
    hist->SetDirectory(0);

    TH1F *zoomedHist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Deposited Energy [keV]; Counts / %d eV",
             Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    zoomedHist->SetDirectory(0);

    TH1F *peakHist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Deposited Energy [keV]; Counts / %d eV",
             Constants::BIN_WIDTH_EV),
        Constants::PEAK_NBINS, Constants::PEAK_XMIN, Constants::PEAK_XMAX);
    peakHist->SetDirectory(0);

    Int_t n_entries = tree->GetEntries();

    for (Int_t j = 0; j < n_entries; j++) {
      tree->GetEntry(j);

      Bool_t in_excluded_region = kFALSE;
      if (nInteractions != 1)
        in_excluded_region = kTRUE;
      if (z < Constants::FILTER_DEPTH_UM)
        in_excluded_region = kTRUE;

      if (!in_excluded_region) {
        for (size_t r = 0;
             r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size(); r++) {
          if (x >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin &&
              x <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax &&
              y >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin &&
              y <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax) {
            in_excluded_region = kTRUE;
            break;
          }
        }
      }

      if (!in_excluded_region) {
        deposited_energy = calibration_function->Eval(energy);
        cal_tree->Fill();
        hist->Fill(deposited_energy);
        zoomedHist->Fill(deposited_energy);
        peakHist->Fill(deposited_energy);
      }
    }

    std::cout << "Calibrated " << input_name << ": " << cal_tree->GetEntries()
              << " / " << n_entries << " events passed filter" << std::endl;

    file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);
    file->Close();
    delete file;
    delete hist;
    delete zoomedHist;
    delete peakHist;
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t use_calibrated = kFALSE;
  Bool_t interactive = kTRUE;

  TF1 *cal_jan13 = BuildCalibration_20260113(use_calibrated, interactive);
  TF1 *cal_jan14 = BuildCalibration_20260114(use_calibrated, interactive);
  TF1 *cal_jan15 = BuildCalibration_20260115(use_calibrated, interactive);
  TF1 *cal_jan16 = BuildCalibration_20260116(use_calibrated, interactive);

  std::vector<TString> datasets_jan13 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113,
      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
      Constants::POSTREACTOR_AM241_20260113,
      Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113,
      Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113};
  PulseHeightToDepositedEnergy(datasets_jan13, cal_jan13, "20260113");

  std::vector<TString> datasets_jan14 = {
      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
      Constants::CUSHIELDSIGNAL_90PERCENT_20260114};
  PulseHeightToDepositedEnergy(datasets_jan14, cal_jan14, "20260114");

  std::vector<TString> datasets_jan15 = {
      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115,
      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
      Constants::POSTREACTOR_AM241_20260115,
      Constants::POSTREACTOR_BA133_20260115,
      Constants::SHUTTERCLOSED_20260115};
  PulseHeightToDepositedEnergy(datasets_jan15, cal_jan15, "20260115");

  std::vector<TString> datasets_jan16 = {
      Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116,
      Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116,
      Constants::POSTREACTOR_AM241_BA133_20260116};
  PulseHeightToDepositedEnergy(datasets_jan16, cal_jan16, "20260116");
}
