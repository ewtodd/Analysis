#include "InitUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void UnfoldData(std::vector<TString> input_filepaths,
                std::vector<TString> output_names, Bool_t reprocess = kFALSE) {
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
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("Timestamp", &timestamp);

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    const Int_t n_boards = 4;
    const Int_t n_channels = 15;
    TTree *trees[n_boards][n_channels];
    UShort_t energies[n_boards][n_channels];
    ULong64_t timestamps[n_boards][n_channels];

    for (Int_t b = 0; b < n_boards; b++) {
      for (Int_t ch = 0; ch < n_channels; ch++) {
        TString tree_name = Form("board%d_ch%d", b, ch);
        trees[b][ch] = new TTree(tree_name, tree_name);
        trees[b][ch]->Branch("energy", &energies[b][ch], "energy/s");
        trees[b][ch]->Branch("timestamp", &timestamps[b][ch], "timestamp/l");
      }
    }

    Int_t n_entries = input_tree->GetEntries();
    for (Int_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);

      if (board < n_boards && channel < n_channels) {
        energies[board][channel] = energy;
        timestamps[board][channel] = timestamp;
        trees[board][channel]->Fill();
      }
    }

    file->Close();

    output_file->cd();
    for (Int_t b = 0; b < n_boards; b++) {
      for (Int_t ch = 0; ch < n_channels; ch++) {
        if (trees[b][ch]->GetEntries() > 0) {
          trees[b][ch]->Write();
        }
        delete trees[b][ch];
      }
    }
    output_file->Close();
  }
}

void EventExtractor() {
  Bool_t reprocess_initial = kTRUE;

  const Int_t n_run_files = 4;
  std::vector<TString> filepaths, output_names;

  TString path_prefix =
      "/home/e-work/LabData/MUSIC/37Cl/ProductionMode/run_37/ROOT/";
  for (Int_t i = 0; i < n_run_files; i++) {
    filepaths.push_back(Form("%sDataR_run_37_%d.root", path_prefix.Data(), i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filepaths[i] << std::endl;
    output_names.push_back(Form("Run37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  UnfoldData(filepaths, output_names, reprocess_initial);
}
