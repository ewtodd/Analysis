
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include <TROOT.h>
#include <vector>

namespace Constants {

const Int_t ADC_MAX = 16384;

// Pulse Height histogram (ADC units)
const Int_t PH_BIN_WIDTH = 20;
const Int_t PH_HIST_XMIN = 0;
const Int_t PH_HIST_XMAX = ADC_MAX;
const Int_t PH_HIST_NBINS = (PH_HIST_XMAX - PH_HIST_XMIN) / PH_BIN_WIDTH;

// Pulse Integral histogram (a.u.)
const Int_t PI_BIN_WIDTH = 800;
const Int_t PI_HIST_XMIN = 0;
const Int_t PI_HIST_XMAX = 400000;
const Int_t PI_HIST_NBINS = (PI_HIST_XMAX - PI_HIST_XMIN) / PI_BIN_WIDTH;

// Light Output histogram (keVee)
const Int_t LO_BIN_WIDTH = 1;
const Float_t LO_HIST_XMIN = 0;
const Float_t LO_HIST_XMAX = 1250;
const Int_t LO_HIST_NBINS = (LO_HIST_XMAX - LO_HIST_XMIN) / LO_BIN_WIDTH;

const Float_t E_LA_33KEV = 33.4;
const Float_t E_AM241_59KEV = 59.5409;
const Float_t E_EU152_122KEV = 121.8;
const Float_t E_EU152_245KEV = 244.7;
const Float_t E_EU152_344KEV = 344.3;

const TString CALIBRATION_AM241 = "calibration_Am241";
const TString CALIBRATION_EU152 = "calibration_Eu152";
const TString BACKGROUND = "background";
const TString IRRADIATION_ONE = "irradiation_one";
const TString IRRADIATION_TWO = "irradiation_two";
const TString IRRADIATION_THREE = "irradiation_three";
const TString IRRADIATION_FOUR = "irradiation_four";

const std::vector<TString> ALL_DATASETS = {
    CALIBRATION_AM241, CALIBRATION_EU152, BACKGROUND,      IRRADIATION_ONE,
    IRRADIATION_TWO,   IRRADIATION_THREE, IRRADIATION_FOUR};

const std::vector<TString> IRRADIATION_DATASETS = {
    IRRADIATION_ONE, IRRADIATION_TWO, IRRADIATION_THREE, IRRADIATION_FOUR};

} // namespace Constants
#endif
