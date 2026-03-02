#include "Constants.hpp"
#include "InitUtils.hpp"
#include <iostream>
#include <vector>

void BinaryToRoot() {
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT);

  Int_t n_files = Constants::ALL_RAW_FILEPATHS.size();

  for (Int_t i = 0; i < n_files; i++) {
    std::cout << "Processing file: " << Constants::ALL_RAW_FILEPATHS[i]
              << std::endl;
    InitUtils::ConvertWavedumpBinToROOT(Constants::ALL_RAW_FILEPATHS[i],
                                        Constants::ALL_OUTPUT_NAMES[i]);
  }
}
