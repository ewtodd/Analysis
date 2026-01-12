#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <vector>

struct TimeShiftResult {
  Long64_t board_0_1_average;
  Long64_t board_0_1_stddev;
  Long64_t board_0_2_average;
  Long64_t board_0_2_stddev;
  Long64_t board_0_3_average;
  Long64_t board_0_3_stddev;
};

std::vector<TimeShiftResult>
CalcAndPlotTimeShifts(std::vector<TString> input_names,
                      Bool_t is_shifted = kFALSE, Bool_t reprocess = kFALSE,
                      Bool_t reprocess_plotting = kFALSE) {
  if (!reprocess)
    return {};

  std::vector<TimeShiftResult> results;
  Int_t n_files = input_names.size();

  for (Int_t i = 0; i < n_files; i++) {
    TString input_name = input_names[i];
    TString input_filepath = "root_files/" + input_name + ".root";
    TFile *input_file = new TFile(input_filepath, "READ");
    const Int_t n_boards = 4;
    TTree *ch5_trees[n_boards];
    ULong64_t ch5_timestamps[n_boards];

    TString branch_name = is_shifted ? "ShiftedTimestamp" : "Timestamp";

    for (Int_t b = 0; b < n_boards; b++) {
      TString tree_name = Form("100HzPulserBoard%d", b);
      ch5_trees[b] = static_cast<TTree *>(input_file->Get(tree_name));
      ch5_trees[b]->SetBranchAddress(branch_name, &ch5_timestamps[b]);
    }

    Int_t n_entries_0 = ch5_trees[0]->GetEntries();
    Int_t n_entries_1 = ch5_trees[1]->GetEntries();
    Int_t n_entries_2 = ch5_trees[2]->GetEntries();
    Int_t n_entries_3 = ch5_trees[3]->GetEntries();

    Int_t n_entries =
        TMath::Min(5000, TMath::Min(TMath::Min(n_entries_0, n_entries_1),
                                    TMath::Min(n_entries_2, n_entries_3)));

    std::vector<Long64_t> timediff_0_1;
    std::vector<Long64_t> timediff_0_2;
    std::vector<Long64_t> timediff_0_3;

    for (Int_t j = 0; j < n_entries; j++) {
      ch5_trees[0]->GetEntry(j);
      ch5_trees[1]->GetEntry(j);
      ch5_trees[2]->GetEntry(j);
      ch5_trees[3]->GetEntry(j);

      timediff_0_1.push_back(ch5_timestamps[0] - ch5_timestamps[1]);
      timediff_0_2.push_back(ch5_timestamps[0] - ch5_timestamps[2]);
      timediff_0_3.push_back(ch5_timestamps[0] - ch5_timestamps[3]);
    }

    TimeShiftResult result;

    result.board_0_1_average =
        TMath::Mean(timediff_0_1.begin(), timediff_0_1.end());
    result.board_0_1_stddev =
        TMath::RMS(timediff_0_1.begin(), timediff_0_1.end());

    result.board_0_2_average =
        TMath::Mean(timediff_0_2.begin(), timediff_0_2.end());
    result.board_0_2_stddev =
        TMath::RMS(timediff_0_2.begin(), timediff_0_2.end());

    result.board_0_3_average =
        TMath::Mean(timediff_0_3.begin(), timediff_0_3.end());
    result.board_0_3_stddev =
        TMath::RMS(timediff_0_3.begin(), timediff_0_3.end());

    results.push_back(result);

    std::cout << input_name << std::endl;
    std::cout << "Board 0, 1 time difference: "
              << result.board_0_1_average * 1e-12 << " s" << std::endl;
    std::cout << "Board 0, 2 time difference: "
              << result.board_0_2_average * 1e-12 << " s" << std::endl;
    std::cout << "Board 0, 3 time difference: "
              << result.board_0_3_average * 1e-12 << " s" << std::endl;
    std::cout << std::endl;

    if (reprocess_plotting) {
      Float_t range_0_1_min =
          (result.board_0_1_average - 5 * result.board_0_1_stddev) * 1e-12;
      Float_t range_0_1_max =
          (result.board_0_1_average + 5 * result.board_0_1_stddev) * 1e-12;
      Float_t range_0_2_min =
          (result.board_0_2_average - 5 * result.board_0_2_stddev) * 1e-12;
      Float_t range_0_2_max =
          (result.board_0_2_average + 5 * result.board_0_2_stddev) * 1e-12;
      Float_t range_0_3_min =
          (result.board_0_3_average - 5 * result.board_0_3_stddev) * 1e-12;
      Float_t range_0_3_max =
          (result.board_0_3_average + 5 * result.board_0_3_stddev) * 1e-12;

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

      for (Int_t j = 0; j < n_entries; j++) {
        hist_0_1->Fill(timediff_0_1[j] * 1e-12);
        hist_0_2->Fill(timediff_0_2[j] * 1e-12);
        hist_0_3->Fill(timediff_0_3[j] * 1e-12);
      }

      TString prefix = is_shifted ? "shiftedtimediff" : "timediff";

      TCanvas *canvas_0_1 =
          new TCanvas(Form("canvas_0_1_%s", input_name.Data()),
                      Form("Time Diff 0-1: %s", input_name.Data()), 1600, 900);

      PlottingUtils::ConfigureCanvas(canvas_0_1, kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(
          hist_0_1, color_0_1,
          Form("%s, Time difference: %.2e", input_name.Data(),
               result.board_0_1_average * 1e-12));
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
               result.board_0_2_average * 1e-12));
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
               result.board_0_3_average * 1e-12));
      PlottingUtils::SaveFigure(
          canvas_0_3, Form("%s_0_3_%s.png", prefix.Data(), input_name.Data()),
          kFALSE);

      delete canvas_0_1;
      delete canvas_0_2;
      delete canvas_0_3;
      delete hist_0_1;
      delete hist_0_2;
      delete hist_0_3;
    }

    input_file->Close();
    delete input_file;
  }

  return results;
}
void ExtractPulserData(std::vector<TString> input_filepaths,
                       std::vector<TString> output_names,
                       Bool_t shifted = kFALSE, Bool_t reprocess = kTRUE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {

    TString input_filepath = input_filepaths[i];
    TFile *file = new TFile(input_filepath, "READ");

    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file for filepath " << input_filepath
                << std::endl;
      continue;
    }

    TTree *input_tree = static_cast<TTree *>(file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree for filepath " << input_filepath
                << std::endl;
      file->Close();
      continue;
    }

    UShort_t board, channel, energy;
    UInt_t flags;
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    if (!shifted)
      input_tree->SetBranchAddress("Timestamp", &timestamp);
    else
      input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    const Int_t n_boards = Constants::N_BOARDS;

    TTree *trees[n_boards];
    ULong64_t timestamps[n_boards];

    for (Int_t b = 0; b < n_boards; b++) {
      TString tree_name = Constants::channelMap.at({b, 5});
      trees[b] = new TTree(tree_name, tree_name);
      if (!shifted)
        trees[b]->Branch("Timestamp", &timestamps[b], "Timestamp/l");
      else
        trees[b]->Branch("ShiftedTimestamp", &timestamps[b],
                         "ShiftedTimestamp/l");
    }

    Int_t n_entries = input_tree->GetEntries();
    for (Int_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      if (channel == 5) {
        timestamps[board] = timestamp;
        trees[board]->Fill();
      }
    }

    file->Close();

    output_file->cd();
    for (Int_t b = 0; b < n_boards; b++) {
      if (trees[b]->GetEntries() > 0) {
        trees[b]->Write();
        delete trees[b];
      }
    }
    output_file->Close();
  }
}

void ApplyTimeShift(std::vector<TString> input_filepaths,
                    std::vector<TimeShiftResult> results, Bool_t reprocess) {
  if (!reprocess)
    return;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filepaths[i];
    TimeShiftResult res = results[i];

    Long64_t mean_0_1, mean_0_2, mean_0_3;
    mean_0_1 = res.board_0_1_average;
    mean_0_2 = res.board_0_2_average;
    mean_0_3 = res.board_0_3_average;

    TFile *input_output_file = new TFile(input_filepath, "UPDATE");

    TTree *input_tree = static_cast<TTree *>(input_output_file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree for filepath " << input_filepath
                << std::endl;
      input_output_file->Close();
      continue;
    }

    UShort_t board, energy;
    UInt_t flags;
    ULong64_t timestamp, shifted_timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Timestamp", &timestamp);
    input_tree->Branch("ShiftedTimestamp", &shifted_timestamp,
                       "ShiftedTimestamp/l");

    Long64_t shifted_timestamp_calc;

    Int_t n_entries = input_tree->GetEntries();
    for (Int_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      switch (board) {
      case 0:
        shifted_timestamp_calc = timestamp;
        break;
      case 1:
        shifted_timestamp_calc = timestamp + mean_0_1;
        break;
      case 2:
        shifted_timestamp_calc = timestamp + mean_0_2;
        break;
      case 3:
        shifted_timestamp_calc = timestamp + mean_0_3;
        break;
      }
      shifted_timestamp = shifted_timestamp_calc;
      input_tree->GetBranch("ShiftedTimestamp")->Fill();
    }

    input_output_file->cd();
    input_tree->Write("", TObject::kOverwrite);
    input_output_file->Close();
  }
}

void TimesortData(std::vector<TString> input_filepaths,
                  std::vector<TString> output_names, Bool_t reprocess = kTRUE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {

    TString input_filepath = input_filepaths[i];

    TFile *input_file = new TFile(input_filepath, "READ");
    TTree *input_tree = (TTree *)input_file->Get("Data_R");

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Total entries: " << n_entries << std::endl;

    ULong64_t shifted_timestamp;
    input_tree->SetBranchAddress("ShiftedTimestamp", &shifted_timestamp);

    std::vector<std::pair<ULong64_t, Long64_t>> timestamp_index;
    timestamp_index.reserve(n_entries);

    std::cout << "Reading timestamps..." << std::endl;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      timestamp_index.push_back({shifted_timestamp, j});

      if (j % 10000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    std::cout << "Sorting..." << std::endl;
    std::sort(timestamp_index.begin(), timestamp_index.end());

    std::cout << "Writing sorted tree..." << std::endl;

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE", "", 0);

    TTree *output_tree = input_tree->CloneTree(0);

    Long64_t *index_array = new Long64_t[n_entries];
    for (Long64_t j = 0; j < n_entries; j++) {
      index_array[j] = timestamp_index[j].second;
    }

    std::cout << "Bulk copying entries..." << std::endl;
    output_tree->CopyEntries(input_tree, -1, "auto", index_array);

    delete[] index_array;

    output_file->cd();
    output_tree->Write();
    output_file->Close();
    input_file->Close();
    std::cout << "Sorted file saved to: " << output_filepath << std::endl;
  }
}

void Timing() {
  Bool_t reprocess_calculation = kTRUE;
  Bool_t reprocess_plotting = kFALSE;
  Bool_t reprocess_sorting = kTRUE;

  InitUtils::SetROOTPreferences();

  std::vector<TString> filepaths, output_names, shifted_output_names,
      sorted_output_names;

  TString path_prefix = "./root_files/";
  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filepaths.push_back(Form("%sDataR_run_37_%d.root", path_prefix.Data(), i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filepaths[i] << std::endl;
    output_names.push_back(Form("TimeShift_Run37_%d", i));
    shifted_output_names.push_back(Form("ShiftedTimeShift_Run37_%d", i));
    sorted_output_names.push_back(Form("Sorted_Run37_%d", i));
  }

  Bool_t is_shifted = kFALSE;
  ExtractPulserData(filepaths, output_names, kFALSE, reprocess_calculation);
  std::vector<TimeShiftResult> results = CalcAndPlotTimeShifts(
      output_names, is_shifted, reprocess_calculation, reprocess_plotting);
  ApplyTimeShift(filepaths, results, kTRUE);
  is_shifted = kTRUE;
  // reprocess_calculation = kFALSE;
  ExtractPulserData(filepaths, shifted_output_names, kTRUE,
                    reprocess_calculation);
  std::vector<TimeShiftResult> check_results =
      CalcAndPlotTimeShifts(shifted_output_names, is_shifted,
                            reprocess_calculation, reprocess_plotting);
  TimesortData(filepaths, sorted_output_names, reprocess_sorting);
}
