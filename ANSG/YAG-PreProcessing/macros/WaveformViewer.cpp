#include "Constants.hpp"
#include "InitUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <iostream>
#include <vector>

void WaveformViewer() {
  InitUtils::SetROOTPreferences(Constants::SAVE_FORMAT);

  WaveformProcessingUtils::ProcessFilesParallel(
      Constants::ALL_FILEPATHS, Constants::ALL_OUTPUT_NAMES,
      Constants::DEFAULT_PROCESSING_CONFIG);
}
