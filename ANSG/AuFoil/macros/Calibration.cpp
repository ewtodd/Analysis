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

TH1D *LoadHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1D *hist = static_cast<TH1D *>(file->Get("zoomedHist"));
  if (!hist) {
    std::cerr << "ERROR: Cannot find zoomedHist in " << input_name << std::endl;
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

  Bool_t use_flat_background = kFALSE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  if (peak_name == "Ba_80.98keV") {
    fitter =
        new FittingUtils(hist, 75, 90, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  } else if (peak_name == "Ba_53keV") {
    fitter =
        new FittingUtils(hist, 45, 57, use_flat_background, use_step,
                         use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);
  }

  if (interactive)
    fitter->SetInteractive();
  FitResult result = fitter->FitSinglePeak(input_name, peak_name);
  delete hist;
  delete fitter;
  return result;
}

FitResult FitPbKAlpha(const TString input_name, const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kFALSE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter = nullptr;

  use_flat_background = kFALSE;
  fitter =
      new FittingUtils(hist, 65, 81, use_flat_background, use_step,
                       use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

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

TF1 *BuildCalibration(const Bool_t interactive) {
  std::cout << "Building Calibration" << std::endl;

  CalibrationData cal_data;
  AddZeroPoint(cal_data);

  FitResult ba_53 =
      FitCalibrationPeak(Constants::YESGOLD, "Ba_53keV", interactive);
  AddCalibrationPoint(cal_data, "Ba_53keV", ba_53.peaks.at(0).mu,
                      ba_53.peaks.at(0).mu_error, E_BA133_53,
                      ba_53.reduced_chi2);

  FitResult bkg_pb = FitPbKAlpha(Constants::YESGOLD, interactive);
  AddCalibrationPoint(cal_data, "Bkg Pb-Ka1", bkg_pb.peaks.at(0).mu,
                      bkg_pb.peaks.at(0).mu_error, E_PB_KA1,
                      bkg_pb.reduced_chi2);
  AddCalibrationPoint(cal_data, "Bkg Pb-Ka2", bkg_pb.peaks.at(1).mu,
                      bkg_pb.peaks.at(1).mu_error, E_PB_KA2, -1);

  FitResult ba_80 =
      FitCalibrationPeak(Constants::YESGOLD, "Ba_80.98keV", interactive);
  AddCalibrationPoint(cal_data, "Ba_80keV", ba_80.peaks.at(0).mu,
                      ba_80.peaks.at(0).mu_error, E_BA133_81,
                      ba_80.reduced_chi2);

  PrintCalibrationSummary(cal_data, "Gold");

  TF1 *cal_func = CreateAndSaveCalibration(
      cal_data.mu, cal_data.calibration_values_keV, cal_data.mu_errors, "Gold");
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

    Double_t deposited_energy = 0;
    TTree *cal_tree = new TTree("calibrated_tree", "Calibrated events");
    cal_tree->SetDirectory(file);
    cal_tree->Branch("depositedEnergykeV", &deposited_energy,
                     "depositedEnergykeV/D");

    TH1D *hist = new TH1D(PlottingUtils::GetRandomName(),
                          Form("; Deposited Energy [keV]; Counts / %d eV",
                               Constants::BIN_WIDTH_EV),
                          Constants::HIST_NBINS, Constants::HIST_XMIN,
                          Constants::HIST_XMAX);
    hist->SetDirectory(0);

    TH1D *zoomedHist =
        new TH1D(PlottingUtils::GetRandomName(),
                 Form("; Deposited Energy [keV]; Counts / %d eV",
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);
    zoomedHist->SetDirectory(0);

    TH1D *peakHist = new TH1D(PlottingUtils::GetRandomName(),
                              Form("; Deposited Energy [keV]; Counts / %d eV",
                                   Constants::BIN_WIDTH_EV),
                              Constants::PEAK_NBINS, Constants::PEAK_XMIN,
                              Constants::PEAK_XMAX);
    peakHist->SetDirectory(0);

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TString treeName = Form("crystal%d_filtered_tree", c);
      TTree *tree = static_cast<TTree *>(file->Get(treeName));
      if (!tree) {
        std::cerr << "ERROR: Could not find " << treeName << " in " << filepath
                  << std::endl;
        continue;
      }

      Double_t energy = 0;
      tree->SetBranchAddress("energykeV", &energy);

      Int_t n_entries = tree->GetEntries();
      for (Int_t j = 0; j < n_entries; j++) {
        tree->GetEntry(j);
        deposited_energy = calibration_function->Eval(energy);
        cal_tree->Fill();
        hist->Fill(deposited_energy);
        zoomedHist->Fill(deposited_energy);
        peakHist->Fill(deposited_energy);
      }
    }

    file->cd();
    cal_tree->Write("calibrated_tree", TObject::kOverwrite);
    hist->Write("calibrated_hist", TObject::kOverwrite);
    zoomedHist->Write("calibrated_zoomedHist", TObject::kOverwrite);
    peakHist->Write("calibrated_peakHist", TObject::kOverwrite);

    delete hist;
    delete zoomedHist;
    delete peakHist;

    std::cout << "Calibrated " << input_name << std::endl;
    file->Close();
    delete file;
  }
}

void Calibration() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  TF1 *cal = BuildCalibration(interactive);

  std::cout << std::endl;
  std::cout << "Calibration Slope: " << std::fixed << std::setprecision(6)
            << cal->GetParameter(1) << " +/- " << cal->GetParError(1)
            << std::endl;

  PulseHeightToDepositedEnergy(Constants::ALL_DATASETS, cal, "Gold");
}
