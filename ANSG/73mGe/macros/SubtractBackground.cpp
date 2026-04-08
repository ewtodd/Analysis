#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1D.h>
#include <TMinuit.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>
#include <cmath>
#include <iomanip>

const Double_t FIT_XMIN = 300;
const Double_t FIT_XMAX = 600;

// Globals for TMinuit FCN
static TH1D *gSignal = nullptr;
static TH1D *gBackground = nullptr;

struct CrystalHistograms {
  TH1D *hist;
  TH1D *zoomedHist;
  TH1D *peakHist;
};

struct Histograms {
  TH1D *hist;
  TH1D *zoomedHist;
  TH1D *peakHist;
  std::vector<CrystalHistograms> crystalHists;
};

struct DatasetPair {
  TString name;
  TString signalFile;
  TString backgroundFile;
};

Double_t ComputeLiveTimePerCrystal(TFile *file, Int_t crystal) {
  TString prefix = Constants::USE_FILTERED ? "Filtered" : "Unfiltered";
  TString paramName = Form("LiveTime_%s_Crystal%d_s", prefix.Data(), crystal);
  TParameter<Double_t> *param =
      (TParameter<Double_t> *)file->Get(paramName);
  if (!param) {
    std::cerr << "WARNING: " << paramName << " not found" << std::endl;
    return 0.0;
  }
  return param->GetVal();
}

Double_t ComputeLiveTime(TFile *file) {
  Double_t totalLiveTime_s = 0;
  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++)
    totalLiveTime_s += ComputeLiveTimePerCrystal(file, c);
  return totalLiveTime_s / Constants::N_CRYSTALS;
}

// Interpolate background histogram at energy E using linear interpolation
Double_t InterpolateHist(TH1D *hist, Double_t E) {
  Int_t bin = hist->FindBin(E);
  Double_t binCenter = hist->GetBinCenter(bin);
  Int_t bin2 = (E >= binCenter) ? bin + 1 : bin - 1;

  if (bin2 < 1 || bin2 > hist->GetNbinsX())
    return hist->GetBinContent(bin);

  Double_t x1 = hist->GetBinCenter(bin);
  Double_t x2 = hist->GetBinCenter(bin2);
  Double_t y1 = hist->GetBinContent(bin);
  Double_t y2 = hist->GetBinContent(bin2);

  return y1 + (y2 - y1) * (E - x1) / (x2 - x1);
}

// Chi2 between signal and A * background(E/gain - delta)
// gain: multiplicative drift, delta: additive shift
Double_t ComputeChi2(TH1D *signal, TH1D *background, Double_t A, Double_t gain,
                     Double_t delta, Int_t &nBins) {
  Int_t binLo = signal->FindBin(FIT_XMIN);
  Int_t binHi = signal->FindBin(FIT_XMAX);

  Double_t chi2 = 0.0;
  nBins = 0;

  for (Int_t bin = binLo; bin <= binHi; bin++) {
    Double_t E = signal->GetBinCenter(bin);
    Double_t s = signal->GetBinContent(bin);
    Double_t b = InterpolateHist(background, E / gain - delta);
    Double_t model = A * b;

    Double_t variance = std::abs(s) + A * A * std::abs(b);
    if (variance <= 0)
      continue;

    Double_t residual = s - model;
    chi2 += residual * residual / variance;
    nBins++;
  }

  return chi2;
}

// TMinuit FCN: par[0] = A, par[1] = gain, par[2] = delta
void MinuitChi2FCN(Int_t & /*npar*/, Double_t * /*grad*/, Double_t &fval,
                   Double_t *par, Int_t /*iflag*/) {
  Int_t nBins = 0;
  fval = ComputeChi2(gSignal, gBackground, par[0], par[1], par[2], nBins);
}

struct FitParams {
  Double_t A;
  Double_t A_error;
  Double_t gain;
  Double_t gain_error;
  Double_t delta;
  Double_t delta_error;
  Double_t chi2ndf;
};

FitParams FitScaleAndShift(TH1D *signal, TH1D *background, Double_t initialA) {
  gSignal = signal;
  gBackground = background;

  TMinuit minuit(3);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(MinuitChi2FCN);

  minuit.DefineParameter(0, "A", initialA, initialA * 0.01, 0, 0);
  minuit.DefineParameter(1, "gain", 1.0, 0.001, 0.9, 1.1);
  minuit.DefineParameter(2, "delta", 0.0, 0.1, -5.0, 5.0);

  minuit.Migrad();

  FitParams result;
  minuit.GetParameter(0, result.A, result.A_error);
  minuit.GetParameter(1, result.gain, result.gain_error);
  minuit.GetParameter(2, result.delta, result.delta_error);

  Int_t nBins = 0;
  Double_t chi2 = ComputeChi2(signal, background, result.A, result.gain,
                              result.delta, nBins);
  Int_t ndf = nBins - 3;
  result.chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;

  gSignal = nullptr;
  gBackground = nullptr;

  return result;
}

Histograms GetRawHistograms(TString filename, Bool_t useCalibrated) {
  Histograms result;

  TFile *file = new TFile("root_files/" + filename + ".root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filename << ".root" << std::endl;
    return {};
  }

  TString prefix = useCalibrated ? "calibrated_" : "";

  for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
    TString histName = Form("%shist_crystal%d", prefix.Data(), c);
    TString zoomedName = Form("%szoomedHist_crystal%d", prefix.Data(), c);
    TString peakName = Form("%speakHist_crystal%d", prefix.Data(), c);

    TH1D *cHist = static_cast<TH1D *>(file->Get(histName));
    TH1D *cZoomed = static_cast<TH1D *>(file->Get(zoomedName));
    TH1D *cPeak = static_cast<TH1D *>(file->Get(peakName));

    if (!cHist || !cZoomed || !cPeak) {
      std::cerr << "WARNING: missing " << prefix << "crystal " << c
                << " histograms for " << filename << std::endl;
      continue;
    }

    CrystalHistograms ch;
    ch.hist = static_cast<TH1D *>(cHist->Clone());
    ch.zoomedHist = static_cast<TH1D *>(cZoomed->Clone());
    ch.peakHist = static_cast<TH1D *>(cPeak->Clone());

    ch.hist->SetDirectory(0);
    ch.zoomedHist->SetDirectory(0);
    ch.peakHist->SetDirectory(0);

    result.crystalHists.push_back(ch);
  }

  file->Close();
  delete file;

  std::cout << filename << ": loaded "
            << (useCalibrated ? "calibrated" : "uncalibrated") << " histograms"
            << std::endl;
  return result;
}

Double_t GetLiveTimeRatio(TString signalFile, TString backgroundFile,
                          Int_t crystal) {
  TFile *sigFile = new TFile("root_files/" + signalFile + ".root", "READ");
  TFile *bkgFile = new TFile("root_files/" + backgroundFile + ".root", "READ");

  Double_t sigLT = ComputeLiveTimePerCrystal(sigFile, crystal);
  Double_t bkgLT = ComputeLiveTimePerCrystal(bkgFile, crystal);

  sigFile->Close();
  bkgFile->Close();
  delete sigFile;
  delete bkgFile;

  if (bkgLT <= 0)
    return 1.0;
  return sigLT / bkgLT;
}

void SubtractBackground() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG);

  Bool_t useCalibrated = kFALSE;

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

  std::cout << "Scale factor fit region: [" << FIT_XMIN << ", " << FIT_XMAX
            << "] keV" << std::endl
            << std::endl;

  std::vector<Histograms> bkgSubtractedResults;

  for (Int_t i = 0; i < (Int_t)datasets.size(); i++) {
    std::cout << "Processing: " << datasets[i].name << std::endl;

    Histograms signal = GetRawHistograms(datasets[i].signalFile, useCalibrated);
    Histograms background =
        GetRawHistograms(datasets[i].backgroundFile, useCalibrated);

    Histograms bkgSubtracted;

    bkgSubtracted.hist = new TH1D(
        PlottingUtils::GetRandomName(),
        Form("hist_%s_BkgSubtracted; Deposited Energy [keV]; Counts / %d eV",
             datasets[i].name.Data(), Constants::BIN_WIDTH_EV),
        Constants::HIST_NBINS, Constants::HIST_XMIN, Constants::HIST_XMAX);
    bkgSubtracted.zoomedHist = new TH1D(
        PlottingUtils::GetRandomName(),
        Form("zoomedHist_%s_BkgSubtracted; Deposited Energy [keV]; Counts / "
             "%d eV",
             datasets[i].name.Data(), Constants::BIN_WIDTH_EV),
        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
        Constants::ZOOMED_XMAX);
    bkgSubtracted.peakHist = new TH1D(
        PlottingUtils::GetRandomName(),
        Form("peakHist_%s_BkgSubtracted; Deposited Energy [keV]; Counts / %d "
             "eV",
             datasets[i].name.Data(), Constants::BIN_WIDTH_EV),
        Constants::PEAK_NBINS, Constants::PEAK_XMIN, Constants::PEAK_XMAX);
    bkgSubtracted.hist->SetDirectory(0);
    bkgSubtracted.zoomedHist->SetDirectory(0);
    bkgSubtracted.peakHist->SetDirectory(0);

    // Sum signal and background across crystals for the fit
    TH1D *sigSum = static_cast<TH1D *>(
        signal.crystalHists[0].zoomedHist->Clone("sigSum_tmp"));
    sigSum->Reset();
    TH1D *bkgSum = static_cast<TH1D *>(
        background.crystalHists[0].zoomedHist->Clone("bkgSum_tmp"));
    bkgSum->Reset();

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      if (c >= (Int_t)signal.crystalHists.size() ||
          c >= (Int_t)background.crystalHists.size())
        continue;
      sigSum->Add(signal.crystalHists[c].zoomedHist);
      bkgSum->Add(background.crystalHists[c].zoomedHist);
    }

    Double_t ltRatio =
        GetLiveTimeRatio(datasets[i].signalFile, datasets[i].backgroundFile, 0);
    FitParams fp = FitScaleAndShift(sigSum, bkgSum, ltRatio);

    std::cout << "  A = " << std::fixed << std::setprecision(6) << fp.A
              << " +/- " << fp.A_error << ", gain = " << std::setprecision(6)
              << fp.gain << " +/- " << fp.gain_error
              << ", delta = " << std::setprecision(4) << fp.delta << " +/- "
              << fp.delta_error << " keV"
              << " (LT ratio = " << std::setprecision(6) << ltRatio
              << ", chi2/ndf = " << std::setprecision(3) << fp.chi2ndf << ")"
              << std::endl;

    delete sigSum;
    delete bkgSum;

    for (Int_t c = 0; c < Constants::N_CRYSTALS; c++) {
      if (c >= (Int_t)signal.crystalHists.size() ||
          c >= (Int_t)background.crystalHists.size()) {
        std::cerr << "WARNING: missing crystal " << c << " histograms for "
                  << datasets[i].name << std::endl;
        continue;
      }

      // Build shifted+scaled background: signal - A * background(E - delta)
      TH1D *subHist = static_cast<TH1D *>(signal.crystalHists[c].hist->Clone(
          Form("hist_%s_crystal%d_BkgSubtracted", datasets[i].name.Data(), c)));
      for (Int_t bin = 1; bin <= subHist->GetNbinsX(); bin++) {
        Double_t E = subHist->GetBinCenter(bin);
        Double_t s = signal.crystalHists[c].hist->GetBinContent(bin);
        Double_t b = InterpolateHist(background.crystalHists[c].hist,
                                     E / fp.gain - fp.delta);
        subHist->SetBinContent(bin, s - fp.A * b);
      }
      subHist->SetDirectory(0);

      TH1D *subZoomed =
          static_cast<TH1D *>(signal.crystalHists[c].zoomedHist->Clone(
              Form("zoomedHist_%s_crystal%d_BkgSubtracted",
                   datasets[i].name.Data(), c)));
      for (Int_t bin = 1; bin <= subZoomed->GetNbinsX(); bin++) {
        Double_t E = subZoomed->GetBinCenter(bin);
        Double_t s = signal.crystalHists[c].zoomedHist->GetBinContent(bin);
        Double_t b = InterpolateHist(background.crystalHists[c].zoomedHist,
                                     E / fp.gain - fp.delta);
        subZoomed->SetBinContent(bin, s - fp.A * b);
      }
      subZoomed->SetDirectory(0);

      TH1D *subPeak =
          static_cast<TH1D *>(signal.crystalHists[c].peakHist->Clone(
              Form("peakHist_%s_crystal%d_BkgSubtracted",
                   datasets[i].name.Data(), c)));
      for (Int_t bin = 1; bin <= subPeak->GetNbinsX(); bin++) {
        Double_t E = subPeak->GetBinCenter(bin);
        Double_t s = signal.crystalHists[c].peakHist->GetBinContent(bin);
        Double_t b = InterpolateHist(background.crystalHists[c].peakHist,
                                     E / fp.gain - fp.delta);
        subPeak->SetBinContent(bin, s - fp.A * b);
      }
      subPeak->SetDirectory(0);

      CrystalHistograms ch = {subHist, subZoomed, subPeak};
      bkgSubtracted.crystalHists.push_back(ch);

      bkgSubtracted.hist->Add(subHist);
      bkgSubtracted.zoomedHist->Add(subZoomed);
      bkgSubtracted.peakHist->Add(subPeak);
    }

    bkgSubtractedResults.push_back(bkgSubtracted);
  }

  TH1D *hist_Combined =
      static_cast<TH1D *>(bkgSubtractedResults[0].hist->Clone("hist_Combined"));
  hist_Combined->Reset();

  TH1D *zoomedHist_Combined = static_cast<TH1D *>(
      bkgSubtractedResults[0].zoomedHist->Clone("zoomedHist_Combined"));
  zoomedHist_Combined->Reset();

  TH1D *peakHist_Combined = static_cast<TH1D *>(
      bkgSubtractedResults[0].peakHist->Clone("peakHist_Combined"));
  peakHist_Combined->Reset();

  for (Int_t i = 0; i < (Int_t)bkgSubtractedResults.size(); i++) {
    hist_Combined->Add(bkgSubtractedResults[i].hist);
    zoomedHist_Combined->Add(bkgSubtractedResults[i].zoomedHist);
    peakHist_Combined->Add(bkgSubtractedResults[i].peakHist);
  }

  hist_Combined->SetTitle(
      Form("Background Subtracted (Scaled); Energy [keV]; Counts / %d eV",
           Constants::BIN_WIDTH_EV));

  zoomedHist_Combined->SetTitle(Form(
      "Background Subtracted (Scaled, Zoomed); Energy [keV]; Counts / %d eV",
      Constants::BIN_WIDTH_EV));

  peakHist_Combined->SetTitle(Form(
      "Background Subtracted (Scaled, Peak Region); Energy [keV]; Counts / %d eV",
      Constants::BIN_WIDTH_EV));

  TCanvas *canvasFull = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->SetFillStyle(0);
  hist_Combined->SetLineWidth(2);
  hist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasFull, "background_subtracted",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TCanvas *canvasZoomed = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->SetFillStyle(0);
  zoomedHist_Combined->SetLineWidth(2);
  zoomedHist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasZoomed, "background_subtracted_zoomed",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TCanvas *canvasPeak = PlottingUtils::GetConfiguredCanvas();
  PlottingUtils::ConfigureHistogram(peakHist_Combined, kP10Violet);
  peakHist_Combined->SetFillStyle(0);
  peakHist_Combined->SetLineWidth(2);
  peakHist_Combined->Draw("HIST");
  PlottingUtils::SaveFigure(canvasPeak, "background_subtracted_peak",
                            "backgroundSubtraction", PlotSaveOptions::kLOG);

  TString outputName = "Combined_BkgSubtracted_Normalized";
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
