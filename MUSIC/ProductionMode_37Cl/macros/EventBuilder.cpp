#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void ResetEvent(UShort_t leftdE[16], UShort_t rightdE[16], UShort_t totaldE[16],
                ULong64_t timestamp_distribution[36], UInt_t all_flags[36],
                Int_t hits[36], UShort_t &cathode, UShort_t &strip0,
                UShort_t &strip17, UShort_t &grid) {

  for (Int_t k = 0; k < 16; k++) {
    leftdE[k] = rightdE[k] = totaldE[k] = 0;
  }
  for (Int_t k = 0; k < 36; k++) {
    timestamp_distribution[k] = 0;
    all_flags[k] = 0;
    hits[k] = 0;
  }
  cathode = strip0 = strip17 = grid = 0;
}

Int_t GetStripNumber(TString map_name) {
  TString num_str = map_name;
  num_str.Remove(0, 1);
  std::cout << num_str.Atoi() << std::endl;
  return num_str.Atoi();
}

void BuildEvents(std::vector<TString> input_filepaths,
                 std::vector<TString> output_names, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  const ULong64_t TIME_WINDOW = Constants::COINCIDENCE_WINDOW;
  std::cout << "Using coincidence window of " << TIME_WINDOW * 1e-6
            << " microseconds..." << std::endl;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {

    TString input_filepath = input_filepaths[i];
    TFile *input_file = new TFile(input_filepath, "READ");

    if (!input_file || input_file->IsZombie()) {
      std::cerr << "Error opening file: " << input_filepath << std::endl;
      continue;
    }

    TTree *input_tree = static_cast<TTree *>(input_file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree from: " << input_filepath << std::endl;
      input_file->Close();
      continue;
    }

    UShort_t board, channel, energy;
    UInt_t flags;
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);
    input_tree->SetBranchAddress("Flags", &flags);

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    UShort_t leftdE[16], rightdE[16], totaldE[16];
    ULong64_t timestamp_distribution[36];
    UInt_t all_flags[36];
    Int_t hits[36];
    UShort_t cathode, strip0, strip17, grid;

    TTree *output_tree = new TTree("event", "MUSIC events");
    output_tree->Branch("LeftdE", leftdE, "LeftdE[16]/s");
    output_tree->Branch("RightdE", rightdE, "RightdE[16]/s");
    output_tree->Branch("TotaldE", totaldE, "TotaldE[16]/s");
    output_tree->Branch("TimestampDistribution", timestamp_distribution,
                        "TimestampDistribution[36]/l");
    output_tree->Branch("AllFlags", all_flags, "AllFlags[36]/i");
    output_tree->Branch("Hits", hits, "Hits[36]/I");
    output_tree->Branch("Cathode", &cathode, "Cathode/s");
    output_tree->Branch("Strip0", &strip0, "Strip0/s");
    output_tree->Branch("Strip17", &strip17, "Strip17/s");
    output_tree->Branch("Grid", &grid, "Grid/s");

    ULong64_t event_time_zero = 0;
    Bool_t in_event = kFALSE;

    ResetEvent(leftdE, rightdE, totaldE, timestamp_distribution, all_flags,
               hits, cathode, strip0, strip17, grid);

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Building events from " << n_entries << " entries..."
              << std::endl;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      TString map_name = Constants::channelMap.at({board, channel});
    }

    output_file->cd();
    output_tree->Write();
    output_file->Close();
    input_file->Close();

    std::cout << "Events saved to: " << output_filepath << std::endl;
  }
}

void EventBuilder() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filepaths, output_names;

  TString path_prefix = "./root_files/";
  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filepaths.push_back(Form("%sDataR_run_37_%d.root", path_prefix.Data(), i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filepaths[i] << std::endl;
    output_names.push_back(Form("Events_Run37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  BuildEvents(filepaths, output_names, reprocess_initial);
}
