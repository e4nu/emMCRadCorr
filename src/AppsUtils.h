#ifndef _APPS_UTILS_H_
#define _APPS_UTILS_H_

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include <iomanip>

namespace e4nu {
  namespace utils {
    double GetCLAS6TargetThickness( double tgt ) ; // Constants set for CLAS6
    std::string GetArg(std::string op, int argc, char ** argv );
    bool ExistArg(std::string op, int argc, char ** argv );
  }
}

#endif



