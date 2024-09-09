#include "PID.h"

Int_t PID::OriginalHits(const string& file, const string& tree)
{
    DisableImplicitMT();
    RDataFrame* dm = new RDataFrame(tree, file);
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "hit_" + outname;

    const vector<string> columns = { "RunNumber", "EventNumber", "EventID_High", "EventID_Low",    // Run and event ID
                                     "CellID", "Hit_Energy", "Hit_X", "Hit_Y", "Hit_Z",    // Hit information
                                     "Acts_RecTrk_No", "Acts_RecTrk_P0", "Acts_TagTrk_No", "Acts_TagTrk_P0", "Missing2_Trk_P",    // Tracker
                                     "ECAL_ClusterSub_N", "ECAL_Cluster_N", "clus_10", "clus_10_tot", "clus_20", "clus_20_tot", "clus_sub_10", "clus_sub_20", "clus_sub_match",    // ECAL: cluster number
                                     "Edep", "Eclus_max", "Eclus_second", "Eclus_max_sec_diff", "Eclus_max_sec_dist",    // ECAL: energy
                                     "Edep_HCAL", "Ecell_max_HCAL", "Edep_SideHCAL", "Ecell_max_SideHCAL" };    // HCAL

    auto fout = dm->Define("ntotal", "(Int_t) ECAL_ECell_XYZ.size()")
    .Define("EventID_High", "(Long_t) EventNumber / 100000")
    .Define("EventID_Low", "(Long_t) EventNumber % 100000")
    .Define("Acts_TagTrk_P0", "(Acts_TagTrk_P.size() > 0) ? TMath::Abs(Acts_TagTrk_P[0]) : 0")
    .Define("Acts_RecTrk_P0", "(Acts_RecTrk_P.size() > 0) ? TMath::Abs(Acts_RecTrk_P[0]) : 0")
    .Define("Missing2_Trk_P", "Acts_TagTrk_P0 - Acts_RecTrk_P0")
    .Define("Edep_HCAL", "(HCAL_E_total.size() > 0) ? HCAL_E_total[0] : 0")
    .Define("Edep_SideHCAL", "(SideHCAL_E_total.size() > 0) ? SideHCAL_E_total[0] : 0")
    .Define("Ecell_max_HCAL", "(HCAL_E_Max_Cell.size() > 0) ? HCAL_E_Max_Cell[0] : 0")
    .Define("Ecell_max_SideHCAL", "(SideHCAL_E_Max_Cell.size() > 0) ? SideHCAL_E_Max_Cell[0] : 0")
    .Define("Hit_X", [] (const vector<Double_t>& ECAL_ECell_XYZ, const Int_t& ntotal)->vector<Double_t>
    {
        vector<Double_t> Hit_X;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            const Int_t layer = i / (nCellsX * nCellsY);
            const Double_t x = (i % nCellsX - nCellsXBias + (0.25 - 0.5 * (layer % 2 == 1)) * staggered_x) * CellWidthX;
            Hit_X.emplace_back(x);
        }
        return Hit_X;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("Hit_Y", [] (const vector<Double_t>& ECAL_ECell_XYZ, const Int_t& ntotal)->vector<Double_t>
    {
        vector<Double_t> Hit_Y;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            Int_t layer = i / (nCellsX * nCellsY);
            const Double_t y = ((i % (nCellsX * nCellsY)) / nCellsX - nCellsYBias + (0.25 - 0.5 * (layer % 2 == 1)) * staggered_y) * CellWidthY;
            Hit_Y.emplace_back(y);
        }
        return Hit_Y;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("Hit_Z", [] (const vector<Double_t>& ECAL_ECell_XYZ, const Int_t& ntotal)->vector<Double_t>
    {
        vector<Double_t> Hit_Z;
        for (Int_t i = 0; i < ntotal; ++i)
        {
            if (ECAL_ECell_XYZ.at(i) <= threshold)
                continue;
            const Double_t z = (i / (nCellsX * nCellsY)) * Thick;
            Hit_Z.emplace_back(z);
        }
        return Hit_Z;
    }, {"ECAL_ECell_XYZ", "ntotal"})
    .Define("CellID", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Y, const vector<Double_t>& Hit_Z)->vector<Int_t>
    {
        vector<Int_t> CellID;
        for (Int_t i = 0; i < Hit_X.size(); ++i)
        {
            const Int_t x = round(Hit_X.at(i) / CellWidthX + nCellsXBias);
            const Int_t y = round(Hit_Y.at(i) / CellWidthY + nCellsYBias);
            const Int_t z = round(Hit_Z.at(i) / Thick);
            const Int_t index = 10000 * z + 100 * x + y;
            CellID.emplace_back(index);
        }
        return CellID;
    }, {"Hit_X", "Hit_Y", "Hit_Z"})
    .Define("Hit_Energy", [] (const vector<Double_t>& ECAL_ECell_XYZ, const Int_t& ntotal)->vector<Double_t>
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
    .Define("Edep", [] (const vector<Double_t>& Hit_Energy)->Double_t
    {
        Double_t Edep = 0.0;
        for (const Double_t& i : Hit_Energy)
            Edep += i;
        return Edep;
    }, {"Hit_Energy"})
    .Define("Eclus_max", "(ECAL_Cluster_N >= 1) ? ECAL_Cluster_E[0] : 0.0")
    .Define("Eclus_second", "(ECAL_Cluster_N >= 2) ? ECAL_Cluster_E[1] : 0.0")
    .Define("Eclus_max_sec_diff", "Eclus_max - Eclus_second")
    .Define("Eclus_max_sec_dist", "(ECAL_Cluster_N >= 2) ? 0.1 * TMath::Sqrt(TMath::Power(ECAL_Cluster_X[0] - ECAL_Cluster_X[1], 2) + TMath::Power(ECAL_Cluster_Y[0] - ECAL_Cluster_Y[1], 2) + TMath::Power(ECAL_Cluster_Z[0] - ECAL_Cluster_Z[1], 2)) : 0.0")    // Convert from mm to cm
    .Define("clus_10", [] (const vector<Double_t>& ECAL_Cluster_E, const Int_t& ECAL_Cluster_N, const Double_t& Eclus_max)->Int_t
    {
        Int_t clus_10 = ECAL_Cluster_N;
        for (Int_t i = 0; i < ECAL_Cluster_N; ++i)
            if (ECAL_Cluster_E.at(i) <= 0.1 * Eclus_max)
                --clus_10;
        return clus_10;
    }, {"ECAL_Cluster_E", "ECAL_Cluster_N", "Eclus_max"})
    .Define("clus_20", [] (const vector<Double_t>& ECAL_Cluster_E, const Int_t& ECAL_Cluster_N, const Double_t& Eclus_max)->Int_t
    {
        Int_t clus_20 = ECAL_Cluster_N;
        for (Int_t i = 0; i < ECAL_Cluster_N; ++i)
            if (ECAL_Cluster_E.at(i) <= 0.2 * Eclus_max)
                --clus_20;
        return clus_20;
    }, {"ECAL_Cluster_E", "ECAL_Cluster_N", "Eclus_max"})
    .Define("clus_sub_10", [] (const vector<Double_t>& ECAL_ClusterSub_E, const Int_t& ECAL_ClusterSub_N)->Int_t
    {
        Int_t clus_sub_10 = ECAL_ClusterSub_N;
        for (Int_t i = 0; i < ECAL_ClusterSub_N; ++i)
            if (ECAL_ClusterSub_E.at(i) <= 0.1 * ECAL_ClusterSub_E.at(0))
                --clus_sub_10;
        return clus_sub_10;
    }, {"ECAL_ClusterSub_E", "ECAL_ClusterSub_N"})
    .Define("clus_sub_20", [] (const vector<Double_t>& ECAL_ClusterSub_E, const Int_t& ECAL_ClusterSub_N)->Int_t
    {
        Int_t clus_sub_20 = ECAL_ClusterSub_N;
        for (Int_t i = 0; i < ECAL_ClusterSub_N; ++i)
            if (ECAL_ClusterSub_E.at(i) <= 0.2 * ECAL_ClusterSub_E.at(0))
                --clus_sub_20;
        return clus_sub_20;
    }, {"ECAL_ClusterSub_E", "ECAL_ClusterSub_N"})
    .Define("clus_sub_match", [] (const vector<Int_t>& ECAL_ClusterSub_matchRecTrk)->Int_t
    {
        Int_t clus_sub_match = 0;
        for (const Int_t& i : ECAL_ClusterSub_matchRecTrk)
            clus_sub_match += (i >= 0);
        return clus_sub_match;
    }, {"ECAL_ClusterSub_matchRecTrk"})
    .Define("clus_10_tot", [] (const vector<Double_t>& ECAL_Cluster_E, const Int_t& ECAL_Cluster_N, const Double_t& Edep)->Int_t
    {
        Int_t clus_10_tot = ECAL_Cluster_N;
        for (Int_t i = 0; i < ECAL_Cluster_N; ++i)
            if (ECAL_Cluster_E.at(i) <= 0.1 * Edep)
                --clus_10_tot;
        return clus_10_tot;
    }, {"ECAL_Cluster_E", "ECAL_Cluster_N", "Edep"})
    .Define("clus_20_tot", [] (const vector<Double_t>& ECAL_Cluster_E, const Int_t& ECAL_Cluster_N, const Double_t& Edep)->Int_t
    {
        Int_t clus_20_tot = ECAL_Cluster_N;
        for (Int_t i = 0; i < ECAL_Cluster_N; ++i)
            if (ECAL_Cluster_E.at(i) <= 0.2 * Edep)
                --clus_20_tot;
        return clus_20_tot;
    }, {"ECAL_Cluster_E", "ECAL_Cluster_N", "Edep"})
    .Snapshot(tree, outname, columns);
    delete dm;

    return 0;
}