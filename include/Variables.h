#ifndef VARIABLES_HH
#define VARIABLES_HH
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/TMVAMultiClassGui.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "Hough.h"
using namespace std;
using namespace ROOT;
using namespace TMath;

class Variables
{
public:
    Variables();
    ~Variables();

    Int_t GenNtuple(const string& file, const string& tree);
};

#endif
