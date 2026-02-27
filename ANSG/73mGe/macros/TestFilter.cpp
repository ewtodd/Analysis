#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void Filter(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree"));

    Float_t energy = 0;
    Float_t x = 0, y = 0, z = 0;

    Int_t nInteractions = 0;

    tree_with_pos->SetBranchAddress("energykeV", &energy);
    tree_with_pos->SetBranchAddress("xum", &x);
    tree_with_pos->SetBranchAddress("yum", &y);
    tree_with_pos->SetBranchAddress("zum", &z);
    if (!Constants::FILTERED)
      tree_with_pos->SetBranchAddress("nInteractions", &nInteractions);

    TH2F *XY = new TH2F(PlottingUtils::GetRandomName(),
                        "; X Position [um]; Y Position [um]", 22, -215, 215, 22,
                        -215, 215);

    TH2F *EZ = new TH2F(PlottingUtils::GetRandomName(),
                        "; Energy [keV]; Interaction Z Position [um]",
                        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                        Constants::ZOOMED_XMAX, 200, 0, 100);

    Int_t n_entries = tree_with_pos->GetEntries();

    TParameter<Double_t> *param =
        (TParameter<Double_t> *)file->Get("N42_RealTime_Total");

    if (!param) {
      std::cerr << "WARNING: Could not find N42_RealTime_Total parameter in "
                << filepath << std::endl;
    }

    TH1F *includedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    TH1F *excludedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    Bool_t in_excluded_region;

    for (Int_t i = 0; i < n_entries; i++) {
      in_excluded_region = kFALSE;

      tree_with_pos->GetEntry(i);

      if (!Constants::FILTERED) {
        if (nInteractions != 1)
          in_excluded_region = kTRUE;
      }

      if (z < Constants::FILTER_DEPTH_UM)
        in_excluded_region = kTRUE;

      if (energy > Constants::ZOOMED_XMIN && energy < Constants::ZOOMED_XMAX) {
        XY->Fill(x, y);
        EZ->Fill(energy, z);

        for (size_t r = 0; r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size();
             r++) {
          if (x >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin &&
              x <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax &&
              y >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin &&
              y <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax) {
            in_excluded_region = kTRUE;
            break;
          }
        }
        if (in_excluded_region) {
          excludedSpectrum->Fill(energy);
        } else
          includedSpectrum->Fill(energy);
      }
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXY = new TCanvas("", "", 1200, 800);
    canvasXY->SetRightMargin(0.2);
    XY->Draw("COLZ");
    PlottingUtils::Configure2DHistogram(XY, canvasXY);
    canvasXY->SetLogz(kFALSE);
    canvasXY->Update();

    for (size_t r = 0; r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size();
         r++) {
      TBox *box = new TBox(Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax);
      box->SetLineColor(kBlack);
      box->SetLineWidth(3);
      box->SetLineStyle(2);
      box->SetFillColorAlpha(kBlack, 0.5);
      box->Draw();
    }

    PlottingUtils::SaveFigure(canvasXY, filename + "_YvsX",
                              PlotSaveOptions::kLINEAR);

    TCanvas *canvasEZ = new TCanvas("", "", 1200, 800);
    canvasEZ->SetRightMargin(0.2);
    EZ->Draw("COLZ");
    PlottingUtils::Configure2DHistogram(EZ, canvasEZ);
    canvasEZ->SetLogz(kFALSE);
    canvasEZ->Update();

    TBox *box = new TBox(Constants::ZOOMED_XMIN, 0, Constants::ZOOMED_XMAX,
                         Constants::FILTER_DEPTH_UM);
    box->SetLineColor(kBlack);
    box->SetLineWidth(3);
    box->SetLineStyle(2);
    box->SetFillColorAlpha(kBlack, 0.5);
    box->Draw();

    PlottingUtils::SaveFigure(canvasEZ, filename + "_ZvsE",
                              PlotSaveOptions::kLINEAR);

    TCanvas *canvasRegions = PlottingUtils::GetConfiguredCanvas(kFALSE);

    Double_t maxY = TMath::Max(includedSpectrum->GetMaximum(),
                               excludedSpectrum->GetMaximum());

    PlottingUtils::ConfigureHistogram(includedSpectrum, kRed);
    includedSpectrum->GetYaxis()->SetRangeUser(0, maxY * 1.1);
    includedSpectrum->SetFillStyle(0);
    includedSpectrum->SetLineWidth(2);
    includedSpectrum->GetYaxis()->SetTitleOffset(1.5);
    includedSpectrum->Draw("HIST");

    PlottingUtils::ConfigureHistogram(excludedSpectrum, kBlack);
    excludedSpectrum->SetFillStyle(0);
    excludedSpectrum->SetLineWidth(2);
    excludedSpectrum->Draw("HIST SAME");

    if (filename == Constants::CDSHIELDSIGNAL_10PERCENT_20260113) {
      TLine *line68 = new TLine(68.752, 0, 68.752, maxY * 1.1);
      line68->SetLineColor(kViolet);
      line68->SetLineWidth(2);
      line68->SetLineStyle(2);
      line68->Draw();

      TLatex *label68 = PlottingUtils::AddText(
          "^{73m}Ge 7/2^{+} #rightarrow 9/2^{+} #gamma ray", 0.47, 0.75);
      label68->SetTextColor(kViolet);
    }

    TLegend *leg = PlottingUtils::AddLegend(0.72, 0.9, 0.75, 0.88);
    leg->AddEntry(includedSpectrum, "Accepted", "l");
    leg->AddEntry(excludedSpectrum, "Rejected", "l");
    leg->Draw();

    canvasRegions->SetLeftMargin(0.2);
    PlottingUtils::SaveFigure(canvasRegions, filename + "_FilterSpectra",
                              PlotSaveOptions::kLINEAR);
    std::cout << "Wrote histograms for " << filename << std::endl;

    XY->Write("XY", TObject::kOverwrite);
    EZ->Write("EZ", TObject::kOverwrite);

    canvasRegions->cd();
    canvasRegions->Write("AcceptedRejectedSpectra", TObject::kOverwrite);

    file->Close();
  }
}

void TestFilter() {
  InitUtils::SetROOTPreferences();

  std::vector<TString> filenames;

  //  filenames.push_back(Constants::PASSIVEBACKGROUND_20260112);
  //  filenames.push_back(Constants::CALIBRATION_20260112);
  //
  //  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_5PERCENT_20260113);
  //  filenames.push_back(Constants::ACTIVEBACKGROUND_TEST_90PERCENT_20260113);
  filenames.push_back(Constants::CDSHIELDSIGNAL_10PERCENT_20260113);
  //  filenames.push_back(Constants::CDSHIELDBACKGROUND_10PERCENT_20260113);
  // filenames.push_back(Constants::CDSHIELDSIGNAL_25PERCENT_20260113);
  //  filenames.push_back(Constants::CDSHIELDBACKGROUND_25PERCENT_20260113);
  //  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260113);
  //  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260113);
  //  filenames.push_back(Constants::POSTREACTOR_AM241_20260113);
  //
  //  filenames.push_back(Constants::CUSHIELDSIGNAL_10PERCENT_20260114);
  //  filenames.push_back(Constants::CUSHIELDBACKGROUND_10PERCENT_20260114);
  //  filenames.push_back(Constants::CUSHIELDSIGNAL_90PERCENT_20260114);
  //
  //  filenames.push_back(Constants::NOSHIELDSIGNAL_5PERCENT_20260115);
  //  filenames.push_back(Constants::NOSHIELDBACKGROUND_5PERCENT_20260115);
  //  filenames.push_back(Constants::POSTREACTOR_AM241_20260115);
  //  filenames.push_back(Constants::POSTREACTOR_BA133_20260115);
  //  filenames.push_back(Constants::SHUTTERCLOSED_20260115);
  //
  //  filenames.push_back(Constants::NOSHIELD_GEONCZT_0_5PERCENT_20260116);
  //  filenames.push_back(Constants::NOSHIELD_ACTIVEBACKGROUND_0_5PERCENT_20260116);
  //  filenames.push_back(
  //      Constants::NOSHIELD_GRAPHITECASTLESIGNAL_10PERCENT_20260116);
  //  filenames.push_back(
  //      Constants::NOSHIELD_GRAPHITECASTLEBACKGROUND_10PERCENT_20260116);
  filenames.push_back(Constants::POSTREACTOR_AM241_BA133_20260116);

  Filter(filenames);
}
