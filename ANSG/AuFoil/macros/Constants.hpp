#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include <TROOT.h>
#include <vector>

struct Region {
  Double_t xmin, xmax, ymin, ymax;
};

namespace Constants {
const Int_t BIN_WIDTH_EV = 75;
const Double_t BIN_WIDTH_KEV = BIN_WIDTH_EV / 1000.0;
const Int_t HIST_XMIN = 0, HIST_XMAX = 1500;
const Int_t ZOOMED_XMIN = 40, ZOOMED_XMAX = 110;
const Int_t PEAK_XMIN = 60, PEAK_XMAX = 80;
const Int_t HIST_NBINS = (HIST_XMAX - HIST_XMIN) / BIN_WIDTH_KEV;
const Int_t ZOOMED_NBINS = (ZOOMED_XMAX - ZOOMED_XMIN) / BIN_WIDTH_KEV;
const Int_t PEAK_NBINS = (PEAK_XMAX - PEAK_XMIN) / BIN_WIDTH_KEV;

const Int_t N_CRYSTALS = 4;
const Double_t TENS_OF_NS_TO_S = 1e-8;

const Bool_t USE_FILTERED = kTRUE;

const Int_t FILTER_DEPTH_MM = 8;
const Double_t PIXEL_ACCEPT_HALFWIDTH_MM = 0.02;
const std::vector<Double_t> PIXEL_CENTERS_X_MM = {
    -20.7690, -18.8370, -16.9470, -15.0570, -13.1670, -11.2350,
    -9.3450,  -7.4550,  -5.5230,  -3.6330,  -1.7430,  1.7430,
    3.6330,   5.5230,   7.4550,   9.3450,   11.2350,  13.1670,
    15.0570,  16.9470,  18.8370,  20.7690};
const std::vector<Double_t> PIXEL_CENTERS_Y_MM = {
    -20.7690, -18.8370, -16.9470, -15.0570, -13.1670, -11.2350,
    -9.3450,  -7.4550,  -5.5230,  -3.6330,  -1.7430,  1.7430,
    3.6330,   5.5230,   7.4550,   9.3450,   11.2350,  13.1670,
    15.0570,  16.9470,  18.8370,  20.7690};

const TString NOGOLD = "NoGold";
const TString YESGOLD = "YesGold";

const std::vector<TString> ALL_DATASETS = {
    NOGOLD,
    YESGOLD,
};

} // namespace Constants
#endif
