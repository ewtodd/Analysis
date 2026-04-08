#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <cmath>
#include <iomanip>
#include <vector>

const Float_t E_AU_KA1 = 68.8037;
const Float_t E_PB_KA1 = 72.8042;
const Float_t E_PB_KA2 = 74.9694;

TH1D *LoadHistogram(const TString input_name) {
  TFile *file = new TFile("root_files/" + input_name + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << input_name << ".root" << std::endl;
    return nullptr;
  }

  TH1D *hist = static_cast<TH1D *>(file->Get("calibrated_zoomedHist"));
  if (!hist) {
    std::cerr << "ERROR: Cannot find calibrated_zoomedHist in " << input_name
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

FitResult FitPbKAlpha(const TString input_name, const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kFALSE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter =
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

FitResult FitAuWithConstrainedPb(const TString input_name,
                                 const FitResult &pb_result,
                                 const Bool_t interactive) {
  TH1D *hist = LoadHistogram(input_name);
  if (!hist)
    return {};

  Bool_t use_flat_background = kTRUE;
  Bool_t use_step = kFALSE;
  Bool_t use_low_exp_tail = kTRUE;
  Bool_t use_low_lin_tail = kTRUE;
  Bool_t use_high_exp_tail = kTRUE;

  FittingUtils *fitter =
      new FittingUtils(hist, 62, 81, use_flat_background, use_step,
                       use_low_exp_tail, use_low_lin_tail, use_high_exp_tail);

  if (interactive)
    fitter->SetInteractive();
  FitResult result =
      fitter->FitTriplePeak(input_name, "Au_KAlpha1", pb_result, E_AU_KA1);
  delete hist;
  delete fitter;
  return result;
}

void Fits() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  FitResult pb_result = FitPbKAlpha(Constants::YESGOLD, interactive);

  FitResult au_result =
      FitAuWithConstrainedPb(Constants::YESGOLD, pb_result, interactive);

  if (!au_result.peaks.empty()) {
    std::cout << std::endl;
    std::cout << "Au K-alpha 1 Result:" << std::endl;
    std::cout << "  mu: " << std::fixed << std::setprecision(4)
              << au_result.peaks.at(0).mu << " +/- "
              << au_result.peaks.at(0).mu_error << " keV" << std::endl;
    std::cout << "  chi2/ndf = " << std::setprecision(3)
              << au_result.reduced_chi2 << std::endl;
  }
}
