#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <vector>

struct TimeShiftResult {
  Long64_t board_0_1_average;
  Long64_t board_0_1_stddev;
  Long64_t board_0_2_average;
  Long64_t board_0_2_stddev;
  Long64_t board_0_3_average;
  Long64_t board_0_3_stddev;
};

std::vector<TimeShiftResult> CalcTimeShifts(std::vector<TString> input_names,
                                            Bool_t check = kFALSE,
                                            Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return {};

  std::vector<TimeShiftResult> results;
  Int_t n_inputs = input_names.size();

  for (Int_t i = 0; i < n_inputs; i++) {
    TString input_name = input_names[i];
    TString input_filepath = "root_files/" + input_name + ".root";
    TFile *input_file = new TFile(input_filepath, "READ");
    const Int_t n_boards = 4;
    TTree *ch5_trees[n_boards];
    ULong64_t ch5_timestamps[n_boards];

    TString branch_name = check ? "shifted_timestamp" : "timestamp";

    for (Int_t b = 0; b < n_boards; b++) {
      TString tree_name = Form("board%d_ch5", b);
      ch5_trees[b] = static_cast<TTree *>(input_file->Get(tree_name));
      ch5_trees[b]->SetBranchAddress(branch_name, &ch5_timestamps[b]);
    }

    TString output_filepath =
        check ? Form("root_files/shiftedtimediffs_%s.root", input_name.Data())
              : Form("root_files/timediffs_%s.root", input_name.Data());
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    Long64_t calc_timediff_0_1, calc_timediff_0_2, calc_timediff_0_3;

    TTree *timediff_tree_0_1 = new TTree("timediff_0_1", "timediff_0_1");
    timediff_tree_0_1->Branch("calc_timediff_0_1", &calc_timediff_0_1,
                              "calc_timediff_0_1/L");

    TTree *timediff_tree_0_2 = new TTree("timediff_0_2", "timediff_0_2");
    timediff_tree_0_2->Branch("calc_timediff_0_2", &calc_timediff_0_2,
                              "calc_timediff_0_2/L");

    TTree *timediff_tree_0_3 = new TTree("timediff_0_3", "timediff_0_3");
    timediff_tree_0_3->Branch("calc_timediff_0_3", &calc_timediff_0_3,
                              "calc_timediff_0_3/L");

    Int_t n_entries_0 = ch5_trees[0]->GetEntries();
    Int_t n_entries_1 = ch5_trees[1]->GetEntries();
    Int_t n_entries_2 = ch5_trees[2]->GetEntries();
    Int_t n_entries_3 = ch5_trees[3]->GetEntries();

    Int_t n_entries =
        TMath::Min(5000, TMath::Min(TMath::Min(n_entries_0, n_entries_1),
                                    TMath::Min(n_entries_2, n_entries_3)));

    for (Int_t j = 0; j < n_entries; j++) {
      ch5_trees[0]->GetEntry(j);
      ch5_trees[1]->GetEntry(j);
      ch5_trees[2]->GetEntry(j);
      ch5_trees[3]->GetEntry(j);

      calc_timediff_0_1 = ch5_timestamps[0] - ch5_timestamps[1];
      timediff_tree_0_1->Fill();

      calc_timediff_0_2 = ch5_timestamps[0] - ch5_timestamps[2];
      timediff_tree_0_2->Fill();

      calc_timediff_0_3 = ch5_timestamps[0] - ch5_timestamps[3];
      timediff_tree_0_3->Fill();
    }

    output_file->cd();
    timediff_tree_0_1->Write();
    timediff_tree_0_2->Write();
    timediff_tree_0_3->Write();

    TimeShiftResult result;

    timediff_tree_0_1->Draw("calc_timediff_0_1", "", "goff");
    result.board_0_1_average = TMath::Mean(timediff_tree_0_1->GetSelectedRows(),
                                           timediff_tree_0_1->GetV1());
    result.board_0_1_stddev = TMath::RMS(timediff_tree_0_1->GetSelectedRows(),
                                         timediff_tree_0_1->GetV1());

    timediff_tree_0_2->Draw("calc_timediff_0_2", "", "goff");
    result.board_0_2_average = TMath::Mean(timediff_tree_0_2->GetSelectedRows(),
                                           timediff_tree_0_2->GetV1());
    result.board_0_2_stddev = TMath::RMS(timediff_tree_0_2->GetSelectedRows(),
                                         timediff_tree_0_2->GetV1());

    timediff_tree_0_3->Draw("calc_timediff_0_3", "", "goff");
    result.board_0_3_average = TMath::Mean(timediff_tree_0_3->GetSelectedRows(),
                                           timediff_tree_0_3->GetV1());
    result.board_0_3_stddev = TMath::RMS(timediff_tree_0_3->GetSelectedRows(),
                                         timediff_tree_0_3->GetV1());

    results.push_back(result);
    std::cout << input_name << std::endl;
    std::cout << "Board 0, 1 time difference: "
              << result.board_0_1_average * 1e-12 << " s" << std::endl;
    std::cout << "Board 0, 2 time difference: "
              << result.board_0_2_average * 1e-12 << " s" << std::endl;
    std::cout << "Board 0, 3 time difference: "
              << result.board_0_3_average * 1e-12 << " s" << std::endl;

    std::cout << std::endl;
    output_file->Close();
    input_file->Close();
  }

  return results;
}

void PlotTimeShifts(std::vector<TString> input_names,
                    std::vector<TimeShiftResult> results, Bool_t check,
                    Bool_t reprocess) {
  if (!reprocess)
    return;
  Int_t n_inputs = input_names.size();

  for (Int_t i = 0; i < n_inputs; i++) {
    TString input_name = input_names[i];
    TString input_filepath =
        check ? Form("root_files/shiftedtimediffs_%s.root", input_name.Data())
              : Form("root_files/timediffs_%s.root", input_name.Data());
    TFile *input_file = new TFile(input_filepath, "UPDATE");

    TTree *timediff_tree_0_1 =
        static_cast<TTree *>(input_file->Get("timediff_0_1"));
    TTree *timediff_tree_0_2 =
        static_cast<TTree *>(input_file->Get("timediff_0_2"));
    TTree *timediff_tree_0_3 =
        static_cast<TTree *>(input_file->Get("timediff_0_3"));

    Long64_t calc_timediff_0_1, calc_timediff_0_2, calc_timediff_0_3;
    timediff_tree_0_1->SetBranchAddress("calc_timediff_0_1",
                                        &calc_timediff_0_1);
    timediff_tree_0_2->SetBranchAddress("calc_timediff_0_2",
                                        &calc_timediff_0_2);
    timediff_tree_0_3->SetBranchAddress("calc_timediff_0_3",
                                        &calc_timediff_0_3);

    TimeShiftResult res = results[i];

    Float_t range_0_1_min =
        (res.board_0_1_average - 5 * res.board_0_1_stddev) * 1e-12;
    Float_t range_0_1_max =
        (res.board_0_1_average + 5 * res.board_0_1_stddev) * 1e-12;
    Float_t range_0_2_min =
        (res.board_0_2_average - 5 * res.board_0_2_stddev) * 1e-12;
    Float_t range_0_2_max =
        (res.board_0_2_average + 5 * res.board_0_2_stddev) * 1e-12;
    Float_t range_0_3_min =
        (res.board_0_3_average - 5 * res.board_0_3_stddev) * 1e-12;
    Float_t range_0_3_max =
        (res.board_0_3_average + 5 * res.board_0_3_stddev) * 1e-12;

    std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
    Int_t color_0_1 = colors[0];
    Int_t color_0_2 = colors[1];
    Int_t color_0_3 = colors[2];

    TH1F *hist_0_1 = new TH1F(Form("hist_0_1_%s", input_name.Data()),
                              ";Time Difference Board 0-1 [s];Counts", 50,
                              range_0_1_min, range_0_1_max);

    TH1F *hist_0_2 = new TH1F(Form("hist_0_2_%s", input_name.Data()),
                              ";Time Difference Board 0-2 [s];Counts", 50,
                              range_0_2_min, range_0_2_max);

    TH1F *hist_0_3 = new TH1F(Form("hist_0_3_%s", input_name.Data()),
                              ";Time Difference Board 0-3 [s];Counts", 50,
                              range_0_3_min, range_0_3_max);

    Int_t n_entries = timediff_tree_0_1->GetEntries();
    for (Int_t j = 0; j < n_entries; j++) {
      timediff_tree_0_1->GetEntry(j);
      hist_0_1->Fill(calc_timediff_0_1 * 1e-12);

      timediff_tree_0_2->GetEntry(j);
      hist_0_2->Fill(calc_timediff_0_2 * 1e-12);

      timediff_tree_0_3->GetEntry(j);
      hist_0_3->Fill(calc_timediff_0_3 * 1e-12);
    }

    TString prefix = check ? "shiftedtimediff" : "timediff";

    TCanvas *canvas_0_1 =
        new TCanvas(Form("canvas_0_1_%s", input_name.Data()),
                    Form("Time Diff 0-1: %s", input_name.Data()), 1600, 900);

    PlottingUtils::ConfigureCanvas(canvas_0_1, kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(
        hist_0_1, color_0_1,
        Form("%s, Time difference: %.2e", input_name.Data(),
             res.board_0_1_average * 1e-12));
    PlottingUtils::SaveFigure(
        canvas_0_1, Form("%s_0_1_%s.png", prefix.Data(), input_name.Data()),
        kFALSE);

    TCanvas *canvas_0_2 =
        new TCanvas(Form("canvas_0_2_%s", input_name.Data()),
                    Form("Time Diff 0-2: %s", input_name.Data()), 1600, 900);

    PlottingUtils::ConfigureCanvas(canvas_0_2, kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(
        hist_0_2, color_0_2,
        Form("%s, Time difference: %.2e", input_name.Data(),
             res.board_0_2_average * 1e-12));
    PlottingUtils::SaveFigure(
        canvas_0_2, Form("%s_0_2_%s.png", prefix.Data(), input_name.Data()),
        kFALSE);

    TCanvas *canvas_0_3 =
        new TCanvas(Form("canvas_0_3_%s", input_name.Data()),
                    Form("Time Diff 0-3: %s", input_name.Data()), 1600, 900);

    PlottingUtils::ConfigureCanvas(canvas_0_3, kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(
        hist_0_3, color_0_3,
        Form("%s, Time difference: %.2e", input_name.Data(),
             res.board_0_3_average * 1e-12));
    PlottingUtils::SaveFigure(
        canvas_0_3, Form("%s_0_3_%s.png", prefix.Data(), input_name.Data()),
        kFALSE);

    input_file->cd();
    timediff_tree_0_1->Write();
    timediff_tree_0_2->Write();
    timediff_tree_0_3->Write();
    hist_0_1->Write("timediff_0_1", TObject::kOverwrite);
    hist_0_2->Write("timediff_0_2", TObject::kOverwrite);
    hist_0_3->Write("timediff_0_3", TObject::kOverwrite);
    input_file->Close();

    delete canvas_0_1;
    delete canvas_0_2;
    delete canvas_0_3;
  }
}

void ApplyTimeShift(std::vector<TString> input_names,
                    std::vector<TimeShiftResult> results, Bool_t reprocess) {
  if (!reprocess)
    return;
  Int_t n_inputs = input_names.size();

  for (Int_t i = 0; i < n_inputs; i++) {
    TString input_name = input_names[i];

    TimeShiftResult res = results[i];

    Long64_t mean_0_1, mean_0_2, mean_0_3;
    mean_0_1 = res.board_0_1_average;
    mean_0_2 = res.board_0_2_average;
    mean_0_3 = res.board_0_3_average;

    TString input_filepath = "root_files/" + input_name + ".root";
    TFile *input_output_file = new TFile(input_filepath, "UPDATE");

    const Int_t n_boards = 4;
    const Int_t n_channels = 15;
    TTree *trees[n_boards][n_channels];
    ULong64_t timestamps[n_boards][n_channels];
    ULong64_t shifted_timestamps[n_boards][n_channels];

    for (Int_t b = 0; b < n_boards; b++) {
      for (Int_t ch = 0; ch < n_channels; ch++) {
        TString tree_name = Form("board%d_ch%d", b, ch);
        trees[b][ch] = static_cast<TTree *>(input_output_file->Get(tree_name));
        if (!trees[b][ch]) {
          continue;
        }
        trees[b][ch]->SetBranchAddress("timestamp", &timestamps[b][ch]);
        trees[b][ch]->Branch("shifted_timestamp", &shifted_timestamps[b][ch],
                             "shifted_timestamp/l");

        TTree *input_tree = trees[b][ch];
        Int_t n_entries = input_tree->GetEntries();
        Long64_t shifted_timestamp_calc;
        for (Int_t j = 0; j < n_entries; j++) {
          input_tree->GetEntry(j);
          switch (b) {
          case 0:
            shifted_timestamp_calc = timestamps[b][ch];
            break;
          case 1:
            shifted_timestamp_calc = timestamps[b][ch] + mean_0_1;
            break;
          case 2:
            shifted_timestamp_calc = timestamps[b][ch] + mean_0_2;
            break;
          case 3:
            shifted_timestamp_calc = timestamps[b][ch] + mean_0_3;
            break;
          }
          shifted_timestamps[b][ch] = shifted_timestamp_calc;
          input_tree->GetBranch("shifted_timestamp")->Fill();
        }
      }
    }

    input_output_file->cd();

    for (Int_t b = 0; b < n_boards; b++) {
      for (Int_t ch = 0; ch < n_channels; ch++) {
        TString tree_name = Form("board%d_ch%d", b, ch);
        if (trees[b][ch])
          trees[b][ch]->Write(tree_name, TObject::kOverwrite);
        delete trees[b][ch];
      }
    }

    input_output_file->Close();
  }
}

void TimeShift() {
  Bool_t reprocess_calculation = kTRUE;
  InitUtils::SetROOTPreferences();

  const Int_t n_run_files = 4;
  std::vector<TString> input_names;

  for (Int_t i = 0; i < n_run_files; i++) {
    input_names.push_back(Form("Run37_%d", i));
  }

  std::vector<TimeShiftResult> results =
      CalcTimeShifts(input_names, kFALSE, reprocess_calculation);
  PlotTimeShifts(input_names, results, kFALSE, kTRUE);
  ApplyTimeShift(input_names, results, kTRUE);
  std::vector<TimeShiftResult> check_results =
      CalcTimeShifts(input_names, kTRUE, kTRUE);
  PlotTimeShifts(input_names, check_results, kTRUE, kTRUE);
}
