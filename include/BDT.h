#ifndef BDT_HH
#define BDT_HH
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/TMVAMultiClassGui.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
using namespace std;
using namespace ROOT;

class BDT
{
public:
    BDT();
    ~BDT();

    void AddTrainSig(const string& file, const string& tree)
    {
        train_sig.insert(pair<TString, TString>(TString(file), TString(tree)));
    }

    void AddTrainBkg(const string& file, const string& tree)
    {
        train_bkg.insert(pair<TString, TString>(TString(file), TString(tree)));
    }

    void AddTestSig(const string& file, const string& tree)
    {
        test_sig.insert(pair<TString, TString>(TString(file), TString(tree)));
    }

    void AddTestBkg(const string& file, const string& tree)
    {
        test_bkg.insert(pair<TString, TString>(TString(file), TString(tree)));
    }

    void AddVar(const string& v, const Char_t& type)
    {
        var.insert(pair<TString, Char_t>(TString(v), type));
    }

    Int_t TrainBDT();
    Int_t BDTNtuple(const string& fname, const string& tname);

    void Clear()
    {
        var.clear();
        train_sig.clear();
        train_bkg.clear();
        test_sig.clear();
        test_bkg.clear();
    }


private:
    map<TString, Char_t> var;
    map<TString, TString> train_sig;
    map<TString, TString> train_bkg;
    map<TString, TString> test_sig;
    map<TString, TString> test_bkg;
};

#endif
