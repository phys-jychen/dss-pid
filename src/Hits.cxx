#include "Hits.h"

const Double_t threshold = 0.0;
const Int_t nCellsX = 21, nCellsY = 21;
const Int_t nCellsXBias = (Int_t) round(0.5 * (nCellsX - 1));
const Int_t nCellsYBias = (Int_t) round(0.5 * (nCellsY - 1));
const Double_t CellX = 2.5, CellY = 2.5, CellZ = 4.0;

const Bool_t staggered = true;

Hits::Hits() {}

Hits::~Hits() {}

Int_t Hits::SaveBranches(const string& file, const string& tree)
{
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "sel_" + outname;

    TFile* f = new TFile((TString) file, "READ");
    TTree* t = f->Get<TTree>((TString) tree);
    t->SetBranchStatus("*", 0);

    vector<TString> remains = { "EventNumber", "ECAL_ECell_XYZ", "RunNumber" };
    for (TString re : remains)
        t->SetBranchStatus(re, 1);

    TFile* fnew = new TFile((TString) outname, "RECREATE");
    TTree* tnew = t->CloneTree();
    tnew->Write(0, TObject::kWriteDelete, 0);
    f->Close();
    fnew->Close();

    return 0;
}

Int_t Hits::OriginalHits(const string& file, const string& tree)
{
    DisableImplicitMT();
    RDataFrame* dm = new RDataFrame(tree, file);
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "hit_" + outname;
    auto fout = dm->Define("ntotal", "(Int_t) ECAL_ECell_XYZ.size()")
    .Define("Hit_X", [] (vector<Double_t> ECAL_ECell_XYZ, Int_t ntotal)
    {
        vector<Double_t> Hit_X;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            Int_t layer = i / (nCellsX * nCellsY);
            Double_t x = (i % nCellsX - nCellsXBias) * CellX;
            if (staggered)
            {
                if (layer % 2 == 0)
                    x += 0.25 * CellX;
                else
                    x -= 0.25 * CellX;
            }
            Hit_X.emplace_back(x);
        }
        return Hit_X;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("Hit_Y", [] (vector<Double_t> ECAL_ECell_XYZ, Int_t ntotal)
    {
        vector<Double_t> Hit_Y;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            Int_t layer = i / (nCellsX * nCellsY);
            Double_t y = ((i % (nCellsX * nCellsY)) / nCellsX - nCellsXBias) * CellY;
            if (staggered)
            {
                if (layer % 2 == 0)
                    y += 0.25 * CellY;
                else
                    y -= 0.25 * CellY;
            }
            Hit_Y.emplace_back(y);
        }
        return Hit_Y;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("Hit_Z", [] (vector<Double_t> ECAL_ECell_XYZ, Int_t ntotal)
    {
        vector<Double_t> Hit_Z;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            Double_t z = (i / (nCellsX * nCellsY)) * CellZ;
            Hit_Z.emplace_back(z);
        }
        return Hit_Z;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("Hit_Theta", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z)
    {
        vector<Double_t> theta = {};
        for (Int_t i = 0; i < Hit_X.size(); i++)
        {
            if (Hit_Z.at(i) == 0)
                theta.emplace_back(PiOver2());
            else
            {
                Double_t rho = Sqrt(Power(Hit_X.at(i), 2) + Power(Hit_Y.at(i), 2));
                Double_t angle = ATan2(rho, Hit_Z.at(i));
                theta.emplace_back(angle);
            }
        }
        return theta;
    }, {"Hit_X", "Hit_Y", "Hit_Z"})
    .Define("Hit_Phi", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y)
    {
        vector<Double_t> phi = {};
        for (Int_t i = 0; i < Hit_X.size(); i++)
        {
            if (Hit_X.at(i) == 0)
            {
                if (Hit_Y.at(i) >= 0)
                    phi.emplace_back(0);
                else
                    phi.emplace_back(Pi());
            }
            else
            {
                Double_t angle = ATan2(Hit_Y.at(i), Hit_X.at(i));
                phi.emplace_back(angle);
            }
        }
        return phi;
    }, {"Hit_X", "Hit_Y"})
    .Define("Hit_Energy", [] (vector<Double_t> ECAL_ECell_XYZ, Int_t ntotal)
    {
        vector<Double_t> Hit_Energy;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            else
                Hit_Energy.emplace_back(ECAL_ECell_XYZ.at(i));
        }
        return Hit_Energy;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Snapshot(tree, outname);
    delete dm;

    TFile* f = new TFile((TString) outname, "READ");
    TTree* t = f->Get<TTree>((TString) tree);
    t->SetBranchStatus("*", 0);
    vector<TString> remains = { "EventNumber", "Hit_Energy", "Hit_Phi", "Hit_Theta", "Hit_X", "Hit_Y", "Hit_Z", "RunNumber" };
    for (TString re : remains)
        t->SetBranchStatus(re, 1);
    TFile* fnew = new TFile((TString) outname, "RECREATE");
    TTree* tnew = t->CloneTree();
    tnew->Write(0, TObject::kWriteDelete, 0);
    f->Close();
    fnew->Close();

    return 0;
}