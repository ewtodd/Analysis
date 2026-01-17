#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <TROOT.h>

namespace Constants {

const Int_t BIN_WIDTH_EV = 175;
const Float_t BIN_WIDTH_KEV = BIN_WIDTH_EV / 1000.0;

const Int_t HIST_XMIN = 0, HIST_XMAX = 1500;
const Int_t ZOOMED_XMIN = 45, ZOOMED_XMAX = 80;
const Int_t HIST_NBINS = (HIST_XMAX - HIST_XMIN) / BIN_WIDTH_KEV;
const Int_t ZOOMED_NBINS = (ZOOMED_XMAX - ZOOMED_XMIN) / BIN_WIDTH_KEV;

const Bool_t NORMALIZE_BY_TIME = kFALSE;
const Bool_t FILTERED = kTRUE;

const Int_t FILTER_WINDOW_X_UM = 100;
const Int_t FILTER_WINDOW_Y_UM = 100;
const Int_t FILTER_DEPTH_UM = 80;

} // namespace Constants

#endif
