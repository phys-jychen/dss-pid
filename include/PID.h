#ifndef PID_HH
#define PID_HH

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectoryFile.h"
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

using std::string;
using std::vector;
using std::map;
using std::unordered_map;
using std::pair;
using std::to_string;
using std::cout;
using std::endl;
using std::ofstream;

using namespace ROOT;
using namespace TMath;

// Energy threshold
const Double_t threshold = 0.0;

// Number of cells
const Int_t nCellsX = 21;
const Int_t nCellsY = 21;
const Int_t nLayer = 11;

// Cell size
const Double_t CellWidthX = 2.5;
const Double_t CellWidthY = 2.5;
const Double_t Thick = 4.0;

const Double_t nCellsXBias = 0.5 * (nCellsX - 1);
const Double_t nCellsYBias = 0.5 * (nCellsY - 1);

// Staggered structure
const Bool_t staggered_x = true;
const Bool_t staggered_y = true;

class PID
{
public:
    PID() = default;

    ~PID() = default;

    static Int_t OriginalHits(const string& file, const string& tree);

    static Int_t GenNtuple(const string& file, const string& tree);

    void AddTrainSig(const string& file, const string& tree, const string& type, const Double_t& weight)
    {
        train_sig.insert(pair<pair<TString, TString>, pair<TString, Double_t>>(pair<TString, TString>(TString(file), TString(tree)), pair(TString(type), weight)));
    }

    void AddTestSig(const string& file, const string& tree, const string& type, const Double_t& weight)
    {
        test_sig.insert(pair<pair<TString, TString>, pair<TString, Double_t>>(pair<TString, TString>(TString(file), TString(tree)), pair(TString(type), weight)));
    }

    void AddVar(const string& v, const Char_t& type)
    {
        var.insert(pair<TString, Char_t>(TString(v), type));
    }

    void AddSpec(const string& v, const string& title)
    {
        spec.insert(pair<TString, TString>(TString(v), TString(title)));
    }

    Int_t TrainBDT();

    static Int_t BDTNtuple(const string& file, const string& tree);

    static Int_t PrintCSV(const string& file);

    void Clear()
    {
        var.clear();
        spec.clear();
        train_sig.clear();
        test_sig.clear();
    }


private:
    static Int_t NewScale(const vector<Double_t>& pos_x, const vector<Double_t>& pos_y, const vector<Double_t>& pos_z, const Int_t& RatioX, const Int_t& RatioY, const Int_t& RatioZ);

    map<TString, Char_t> var;
    map<TString, TString> spec;
    map<pair<TString, TString>, pair<TString, Double_t>> train_sig;
    map<pair<TString, TString>, pair<TString, Double_t>> test_sig;
};

#endif
