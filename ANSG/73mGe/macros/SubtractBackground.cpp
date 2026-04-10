#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "InitUtils.hpp"
#include "InteractiveSubtractionEditor.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <iomanip>

const Int_t SNIP_MIN = 10;
const Int_t SNIP_MAX = 60;
const Float_t RESIDUAL_REGION_LOW = 72.0;
const Float_t RESIDUAL_REGION_HIGH = 90.0;

// Estimate baseline from full-range histogram, subtract matching region
// from target. full_hist and target can be the same pointer.
void SubtractBaseline(TH1F *full_hist, TH1F *target, Int_t snip_iterations) {
  TSpectrum spec;
  TH1 *baseline = spec.Background(full_hist, snip_iterations);
  if (!baseline)
    return;

  for (Int_t bin = 1; bin <= target->GetNbinsX(); bin++) {
    Double_t x = target->GetBinCenter(bin);
    Int_t bl_bin = baseline->FindBin(x);
    Double_t bl_val = baseline->GetBinContent(bl_bin);
    target->SetBinContent(bin, target->GetBinContent(bin) - bl_val);
  }
}

struct CrystalHistograms {
  TH1F *hist;
  TH1F *zoomedHist;
  TH1F *peakHist;
};

struct Histograms {
  TH1F *hist;
  TH1F *zoomedHist;
  TH1F *peakHist;
  std::vector<CrystalHistograms> crystalHists;
};

struct DatasetPair {
  TString name;
  TString signalFile;
  TString backgroundFile;
};

Float_t ComputeLiveTimePerCrystal(TFile *file, Int_t crystal) {
  TString prefix = Constants::USE_FILTERED ? "Filtered" : "Unfiltered";
  TString paramName = Form("LiveTime_%s_Crystal%d_s", prefix.Data(), crystal);
  TParameter<Float_t> *param =
      static_cast<TParameter<Float_t> *>(file->Get(paramName));
  if (!param) {
    std::cerr << "WARNING: " << paramName << " not found" << std::endl;
    return 0.0;
  }
  return param->GetVal();
}

void SubtractBackground() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t interactive = kTRUE;

  std::vector<DatasetPair> datasets;

  datasets.push_back({"CdShield_10_01132026",
                      Constants::CDSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CdShield_25_01132026",
                      Constants::CDSHIELDSIGNAL_25PERCENT_20260113,
                      Constants::CDSHIELDBACKGROUND_25PERCENT_20260113});

  datasets.push_back({"CuShield_10_01132026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260113,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260113});

  datasets.push_back({"CuShield_10_01142026",
                      Constants::CUSHIELDSIGNAL_10PERCENT_20260114,
                      Constants::CUSHIELDBACKGROUND_10PERCENT_20260114});

  datasets.push_back({"NoShield_5_01152026",
                      Constants::NOSHIELDSIGNAL_5PERCENT_20260115,
                      Constants::NOSHIELDBACKGROUND_5PERCENT_20260115});

  datasets.push_back(
      {"NoShield_GraphiteCastle_10_01162026",
       Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116,
       Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116});

  std::vector<Histograms> bkgSubtractedResults;

  for (Int_t i = 0; i < (Int_t)datasets.size(); i++) {
    std::cout << "Processing: " << datasets[i].name << std::endl;

    TFile *sigFile =
        new TFile("root_files/" + datasets[i].signalFile + ".root", "READ");
    TFile *bkgFile =
        new TFile("root_files/" + datasets[i].backgroundFile + ".root", "READ");

    if (!sigFile || sigFile->IsZombie() || !bkgFile || bkgFile->IsZombie()) {
      std::cerr << "ERROR: Cannot open files for " << datasets[i].name
                << std::endl;
      if (sigFile)
        delete sigFile;
      if (bkgFile)
        delete bkgFile;
      continue;
    }

    // Load per-crystal live times and compute LT ratios
    Float_t sigLTs[Constants::N_CRYSTALS];
    Float_t ltRatios[Constants::N_CRYSTALS];

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      sigLTs[c] = ComputeLiveTimePerCrystal(sigFile, c);
      Float_t bkgLT = ComputeLiveTimePerCrystal(bkgFile, c);
      ltRatios[c] = (bkgLT > 0) ? sigLTs[c] / bkgLT : 1.0;

      std::cout << "  Crystal " << c << ": LT ratio = " << std::fixed
                << std::setprecision(4) << ltRatios[c]
                << " (sig = " << sigLTs[c] << " s, bkg = " << bkgLT << " s)"
                << std::endl;
    }

    // Scan SNIP iterations: for each, find optimal correction factor in
    // 72-90 keV via least-squares, then compare minimum achievable residuals
    Int_t best_snip = (SNIP_MIN + SNIP_MAX) / 2;
    Double_t best_residual = 1e30;
    Double_t best_correction = 1.0;

    std::cout << "  Scanning SNIP iterations (" << SNIP_MIN << "-" << SNIP_MAX
              << ")..." << std::endl;

    for (Int_t snip = SNIP_MIN; snip <= SNIP_MAX; snip++) {
      // Build combined baseline-subtracted signal and background zoomed hists
      TH1F *sig_combined =
          new TH1F(PlottingUtils::GetRandomName(), "", Constants::ZOOMED_NBINS,
                   Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);
      sig_combined->SetDirectory(0);

      TH1F *bkg_combined =
          new TH1F(PlottingUtils::GetRandomName(), "", Constants::ZOOMED_NBINS,
                   Constants::ZOOMED_XMIN, Constants::ZOOMED_XMAX);
      bkg_combined->SetDirectory(0);

      for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
        TString histName = Form("calibrated_hist_crystal%d", c);
        TString zoomedName = Form("calibrated_zoomedHist_crystal%d", c);
        TH1F *sFull = static_cast<TH1F *>(sigFile->Get(histName));
        TH1F *sZoomed = static_cast<TH1F *>(sigFile->Get(zoomedName));
        TH1F *bFull = static_cast<TH1F *>(bkgFile->Get(histName));
        TH1F *bZoomed = static_cast<TH1F *>(bkgFile->Get(zoomedName));

        if (!sFull || !sZoomed || !bFull || !bZoomed)
          continue;

        TH1F *sbl =
            static_cast<TH1F *>(sZoomed->Clone(PlottingUtils::GetRandomName()));
        sbl->SetDirectory(0);
        SubtractBaseline(sFull, sbl, snip);
        sbl->Scale(1.0 / sigLTs[c]);
        sig_combined->Add(sbl);
        delete sbl;

        TH1F *bbl =
            static_cast<TH1F *>(bZoomed->Clone(PlottingUtils::GetRandomName()));
        bbl->SetDirectory(0);
        SubtractBaseline(bFull, bbl, snip);
        bbl->Scale(ltRatios[c] / sigLTs[c]);
        bkg_combined->Add(bbl);
        delete bbl;
      }

      // Optimal correction: minimize sum((s_i - scale * b_i)^2) in 72-90
      // => scale = sum(s_i * b_i) / sum(b_i^2)
      Double_t sb = 0, bb = 0;
      for (Int_t bin = 1; bin <= sig_combined->GetNbinsX(); bin++) {
        Double_t x = sig_combined->GetBinCenter(bin);
        if (x >= RESIDUAL_REGION_LOW && x <= RESIDUAL_REGION_HIGH) {
          Double_t s = sig_combined->GetBinContent(bin);
          Double_t b = bkg_combined->GetBinContent(bin);
          sb += s * b;
          bb += b * b;
        }
      }

      Double_t opt_scale = (bb > 0) ? sb / bb : 1.0;

      // Minimum residual with this optimal scale
      Double_t residual = 0;
      for (Int_t bin = 1; bin <= sig_combined->GetNbinsX(); bin++) {
        Double_t x = sig_combined->GetBinCenter(bin);
        if (x >= RESIDUAL_REGION_LOW && x <= RESIDUAL_REGION_HIGH) {
          Double_t s = sig_combined->GetBinContent(bin);
          Double_t b = bkg_combined->GetBinContent(bin);
          Double_t r = s - opt_scale * b;
          residual += r * r;
        }
      }

      if (residual < best_residual) {
        best_residual = residual;
        best_snip = snip;
        best_correction = opt_scale;
      }

      delete sig_combined;
      delete bkg_combined;
    }

    std::cout << "  Best SNIP iterations: " << best_snip
              << " (min residual^2 = " << std::scientific
              << std::setprecision(3) << best_residual
              << ", optimal scale = " << std::fixed << std::setprecision(4)
              << best_correction << ")" << std::endl;

    // Build combined histograms for interactive editor using best SNIP
    TH1F *sigZoomedCombined = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    sigZoomedCombined->SetDirectory(0);

    TH1F *bkgZoomedNormed = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    bkgZoomedNormed->SetDirectory(0);

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TString histName = Form("calibrated_hist_crystal%d", c);
      TString zoomedName = Form("calibrated_zoomedHist_crystal%d", c);
      TH1F *sigFull = static_cast<TH1F *>(sigFile->Get(histName));
      TH1F *sigZoomed = static_cast<TH1F *>(sigFile->Get(zoomedName));
      TH1F *bkgFull = static_cast<TH1F *>(bkgFile->Get(histName));
      TH1F *bkgZoomed = static_cast<TH1F *>(bkgFile->Get(zoomedName));

      if (sigFull && sigZoomed && bkgFull && bkgZoomed) {
        TH1F *sigZoomedBL = static_cast<TH1F *>(
            sigZoomed->Clone(PlottingUtils::GetRandomName()));
        sigZoomedBL->SetDirectory(0);
        SubtractBaseline(sigFull, sigZoomedBL, best_snip);

        TH1F *bkgZoomedBL = static_cast<TH1F *>(
            bkgZoomed->Clone(PlottingUtils::GetRandomName()));
        bkgZoomedBL->SetDirectory(0);
        SubtractBaseline(bkgFull, bkgZoomedBL, best_snip);

        sigZoomedCombined->Add(sigZoomedBL);
        bkgZoomedNormed->Add(bkgZoomedBL, ltRatios[c]);

        delete sigZoomedBL;
        delete bkgZoomedBL;
      }
    }

    // Interactive editor: signal vs correction * LT-normalized background
    Double_t correction = best_correction;
    if (interactive) {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      Double_t editor_result = LaunchSubtractionEditor(
          sigZoomedCombined, bkgZoomedNormed, best_correction, 72.0, 90.0,
          datasets[i].name);
      gROOT->SetBatch(was_batch);
      if (editor_result > 0)
        correction = editor_result;
      std::cout << "  Correction factor: " << std::fixed << std::setprecision(4)
                << correction << std::endl;
    }

    delete sigZoomedCombined;
    delete bkgZoomedNormed;

    // Final pass: apply subtraction per crystal using best SNIP
    Histograms bkgSubtracted;

    bkgSubtracted.hist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
    bkgSubtracted.zoomedHist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    bkgSubtracted.peakHist = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        Constants::PEAK_NBINS, Constants::PEAK_XMIN, Constants::PEAK_XMAX);
    bkgSubtracted.hist->SetDirectory(0);
    bkgSubtracted.zoomedHist->SetDirectory(0);
    bkgSubtracted.peakHist->SetDirectory(0);

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      TString histName = Form("calibrated_hist_crystal%d", c);
      TString zoomedName = Form("calibrated_zoomedHist_crystal%d", c);
      TString peakName = Form("calibrated_peakHist_crystal%d", c);

      TH1F *sigHist = static_cast<TH1F *>(sigFile->Get(histName));
      TH1F *sigZoomed = static_cast<TH1F *>(sigFile->Get(zoomedName));
      TH1F *sigPeak = static_cast<TH1F *>(sigFile->Get(peakName));
      TH1F *bkgHist = static_cast<TH1F *>(bkgFile->Get(histName));
      TH1F *bkgZoomed = static_cast<TH1F *>(bkgFile->Get(zoomedName));
      TH1F *bkgPeak = static_cast<TH1F *>(bkgFile->Get(peakName));

      if (!sigHist || !sigZoomed || !sigPeak || !bkgHist || !bkgZoomed ||
          !bkgPeak) {
        std::cerr << "WARNING: missing crystal " << c << " histograms for "
                  << datasets[i].name << std::endl;
        continue;
      }

      Float_t scale = correction * ltRatios[c];

      // Clone all signal histograms and subtract baseline estimated from full
      TH1F *subHist = static_cast<TH1F *>(sigHist->Clone(
          Form("hist_%s_crystal%d_BkgSubtracted", datasets[i].name.Data(), c)));
      subHist->SetDirectory(0);
      SubtractBaseline(sigHist, subHist, best_snip);

      TH1F *subZoomed = static_cast<TH1F *>(
          sigZoomed->Clone(Form("zoomedHist_%s_crystal%d_BkgSubtracted",
                                datasets[i].name.Data(), c)));
      subZoomed->SetDirectory(0);
      SubtractBaseline(sigHist, subZoomed, best_snip);

      TH1F *subPeak = static_cast<TH1F *>(sigPeak->Clone(Form(
          "peakHist_%s_crystal%d_BkgSubtracted", datasets[i].name.Data(), c)));
      subPeak->SetDirectory(0);
      SubtractBaseline(sigHist, subPeak, best_snip);

      // Clone all background histograms and subtract baseline estimated from
      // full
      TH1F *bkgHistBL =
          static_cast<TH1F *>(bkgHist->Clone(PlottingUtils::GetRandomName()));
      bkgHistBL->SetDirectory(0);
      SubtractBaseline(bkgHist, bkgHistBL, best_snip);

      TH1F *bkgZoomedBL =
          static_cast<TH1F *>(bkgZoomed->Clone(PlottingUtils::GetRandomName()));
      bkgZoomedBL->SetDirectory(0);
      SubtractBaseline(bkgHist, bkgZoomedBL, best_snip);

      TH1F *bkgPeakBL =
          static_cast<TH1F *>(bkgPeak->Clone(PlottingUtils::GetRandomName()));
      bkgPeakBL->SetDirectory(0);
      SubtractBaseline(bkgHist, bkgPeakBL, best_snip);

      // Subtract scaled background and normalize to rate
      subHist->Add(bkgHistBL, -scale);
      subHist->Scale(1.0 / sigLTs[c]);

      subZoomed->Add(bkgZoomedBL, -scale);
      subZoomed->Scale(1.0 / sigLTs[c]);

      subPeak->Add(bkgPeakBL, -scale);
      subPeak->Scale(1.0 / sigLTs[c]);

      delete bkgHistBL;
      delete bkgZoomedBL;
      delete bkgPeakBL;

      CrystalHistograms ch = {subHist, subZoomed, subPeak};
      bkgSubtracted.crystalHists.push_back(ch);

      bkgSubtracted.hist->Add(subHist);
      bkgSubtracted.zoomedHist->Add(subZoomed);
      bkgSubtracted.peakHist->Add(subPeak);
    }

    sigFile->Close();
    bkgFile->Close();
    delete sigFile;
    delete bkgFile;

    TCanvas *canvasPair = PlottingUtils::GetConfiguredCanvas();
    PlottingUtils::ConfigureHistogram(bkgSubtracted.peakHist, kP10Violet);
    bkgSubtracted.peakHist->Sumw2(0);
    bkgSubtracted.peakHist->SetMarkerStyle(20);
    bkgSubtracted.peakHist->SetMarkerSize(0.65);
    bkgSubtracted.peakHist->Draw("P");
    PlottingUtils::SaveFigure(canvasPair, "peak_" + datasets[i].name,
                              "backgroundSubtraction/pairs",
                              PlotSaveOptions::kLINEAR);

    bkgSubtractedResults.push_back(bkgSubtracted);
  }

  TH1F *hist_Combined =
      static_cast<TH1F *>(bkgSubtractedResults[0].hist->Clone("hist_Combined"));
  hist_Combined->Reset();

  TH1F *zoomedHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].zoomedHist->Clone("zoomedHist_Combined"));
  zoomedHist_Combined->Reset();

  TH1F *peakHist_Combined = static_cast<TH1F *>(
      bkgSubtractedResults[0].peakHist->Clone("peakHist_Combined"));
  peakHist_Combined->Reset();

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    hist_Combined->Add(bkgSubtractedResults[i].hist);
    zoomedHist_Combined->Add(bkgSubtractedResults[i].zoomedHist);
    peakHist_Combined->Add(bkgSubtractedResults[i].peakHist);
  }

  hist_Combined->SetTitle(
      Form("Background Subtracted; Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));

  zoomedHist_Combined->SetTitle(
      Form("Background Subtracted (Zoomed); Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));

  peakHist_Combined->SetTitle(
      Form("Background Subtracted (Peak Region); Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));

  TCanvas *canvasFull = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->SetMarkerStyle(20);
  hist_Combined->SetMarkerSize(0.65);
  hist_Combined->SetMarkerColor(kP10Violet);
  hist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasFull, "background_subtracted",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  TCanvas *canvasZoomed = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->SetMarkerStyle(20);
  zoomedHist_Combined->SetMarkerSize(0.65);
  zoomedHist_Combined->SetMarkerColor(kP10Violet);
  zoomedHist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasZoomed, "background_subtracted_zoomed",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  TCanvas *canvasPeak = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(peakHist_Combined, kP10Violet);
  peakHist_Combined->Sumw2(0);
  peakHist_Combined->SetMarkerStyle(20);
  peakHist_Combined->SetMarkerSize(0.65);
  peakHist_Combined->SetMarkerColor(kP10Violet);
  peakHist_Combined->Draw("PE");
  PlottingUtils::SaveFigure(canvasPeak, "background_subtracted_peak",
                            "backgroundSubtraction", PlotSaveOptions::kLINEAR);

  // Fit Ge-73m 68.75 keV peak in the combined spectrum (rebinned to 300 eV)
  std::cout << std::endl;
  std::cout
      << "  Fitting Ge-73m peak in combined background-subtracted spectrum:"
      << std::endl;

  FittingUtils *ge_fitter = new FittingUtils(
      zoomedHist_Combined, 64.0, 77.0, kTRUE, kFALSE, kTRUE, kTRUE, kTRUE);
  if (interactive)
    ge_fitter->SetInteractive();
  FitResult ge_fit = ge_fitter->FitSinglePeak("Combined", "Ge73m_68keV");
  if (ge_fit.valid) {
    std::cout << "    mu = " << std::fixed << std::setprecision(4)
              << ge_fit.peaks[0].mu << " +/- " << ge_fit.peaks[0].mu_error
              << " keV, sigma = " << ge_fit.peaks[0].sigma << " +/- "
              << ge_fit.peaks[0].sigma_error
              << " keV, chi2/ndf = " << std::setprecision(3)
              << ge_fit.reduced_chi2 << std::endl;
  } else {
    std::cout << "    FIT FAILED" << std::endl;
  }
  delete ge_fitter;
  std::cout << std::endl;

  TString outputName = "Combined_BkgSubtracted";
  TFile *outputFile =
      new TFile("root_files/" + outputName + ".root", "RECREATE");

  hist_Combined->Write("histCombined", TObject::kOverwrite);
  zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
  peakHist_Combined->Write("peakHistCombined", TObject::kOverwrite);

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    bkgSubtractedResults[i].hist->Write();
    bkgSubtractedResults[i].zoomedHist->Write();
    bkgSubtractedResults[i].peakHist->Write();

    for (Int_t c = 0; c < (Int_t)bkgSubtractedResults[i].crystalHists.size();
         c++) {
      bkgSubtractedResults[i].crystalHists[c].hist->Write();
      bkgSubtractedResults[i].crystalHists[c].zoomedHist->Write();
      bkgSubtractedResults[i].crystalHists[c].peakHist->Write();
    }
  }

  outputFile->Close();
  delete outputFile;

  std::cout << "Background subtraction and combination complete!" << std::endl;
  std::cout << "Output saved to: root_files/" << outputName << ".root"
            << std::endl;
}
