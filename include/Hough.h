#ifndef HOUGH_HH
#define HOUGH_HH
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <ROOT/RDataFrame.hxx>
using namespace std;

class Hough
{
public:
    Hough(const vector<Double_t>& x, const vector<Double_t>& y, const vector<Double_t>& z, const vector<Double_t>& e);
    ~Hough();

    Int_t GetNtrack() {return ntrack;};
    vector<Double_t> GetHclX() {return hcx;};
    vector<Double_t> GetHclY() {return hcy;};
    vector<Double_t> GetHclZ() {return hcz;};
    vector<Double_t> GetHclE() {return hce;};


private:
    Int_t ntrack = 0;

    typedef struct hit
    {
        Double_t x = 0.0;
        Double_t y = 0.0;
        Double_t z = 0.0;
        Double_t e = 0.0;
    }hit;

    typedef struct hcl
    {
        Double_t x = 0.0;
        Double_t y = 0.0;
        Double_t z = 0.0;
        Double_t e = 0.0;
    }hcl;

    TH2D* hht;

    unordered_map<Int_t,vector<hit*>> layer_hits;
    vector<hcl*> Hough_Cluster;
    vector<vector<hit*>> MergeAdjacent(vector<hit*> hits);

    vector<Double_t> hcx;
    vector<Double_t> hcy;
    vector<Double_t> hcz;
    vector<Double_t> hce;
};

#endif
