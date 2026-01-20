#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <TROOT.h>

struct Region {
  Double_t xmin, xmax, ymin, ymax;
};

namespace Constants {

const Int_t BIN_WIDTH_EV = 275;
const Float_t BIN_WIDTH_KEV = BIN_WIDTH_EV / 1000.0;

const Int_t HIST_XMIN = 0, HIST_XMAX = 1500;
const Int_t ZOOMED_XMIN = 60, ZOOMED_XMAX = 80;
const Int_t HIST_NBINS = (HIST_XMAX - HIST_XMIN) / BIN_WIDTH_KEV;
const Int_t ZOOMED_NBINS = (ZOOMED_XMAX - ZOOMED_XMIN) / BIN_WIDTH_KEV;

const Bool_t USE_REAL_TIME = kFALSE;
const Bool_t NORMALIZE_BY_TIME = kTRUE;
const Bool_t FILTERED = kTRUE;

const Int_t FILTER_DEPTH_UM = 80;
const std::vector<Region> FILTER_REGIONS_EXCLUDE_XY_UM = {
    {-215, 215, -215, -205},
    {-215, 215, 205, 215},
    {-215, -205, -205, 205},
    {205, 215, -205, 205}};
} // namespace Constants

#endif
