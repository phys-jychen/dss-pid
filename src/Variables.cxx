#include "PID.h"
#include "Fit3D.h"

Int_t PID::NewScale(const vector<Double_t>& hit_x, const vector<Double_t>& hit_y, const vector<Double_t>& hit_z, const Int_t& ratio_x, const Int_t& ratio_y, const Int_t& ratio_z)
{
    const Int_t nhits = hit_x.size();
    assert(nhits == hit_y.size());
    assert(nhits == hit_z.size());

    unordered_map<Int_t, Int_t> ID_hit_map;

    for (Int_t i = 0; i < nhits; ++i)
    {
        const Double_t x = hit_x.at(i);
        const Double_t y = hit_y.at(i);
        const Int_t z = round(hit_z.at(i) / Thick);

        const Int_t supercellID_x = (Int_t) round(x / CellWidthX + nCellsXBias - (0.25 - 0.5 * (z % 2 == 1)) * staggered_x) / ratio_x;
        const Int_t supercellID_y = (Int_t) round(y / CellWidthY + nCellsYBias - (0.25 - 0.5 * (z % 2 == 1)) * staggered_y) / ratio_y;
        const Int_t supercellID_z = z / ratio_z;
        const Int_t CellID_new = 10000 * supercellID_z + 100 * supercellID_y + supercellID_x;

        ++ID_hit_map[CellID_new];
    }

    const Int_t nhits_new = ID_hit_map.size();
    ID_hit_map.clear();
    return nhits_new;
}

Int_t PID::GenNtuple(const string& file, const string& tree)
{
//    EnableImplicitMT();
    DisableImplicitMT();
    RDataFrame* dm = new RDataFrame(tree, file);
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "rec_" + outname;
    auto fout = dm->Define("nhits", "(Int_t) Hit_X.size()")
    // Transfer z position to the layer number
    .Define("layer", [] (const vector<Double_t>& Hit_Z, const Int_t& nhits)->vector<Int_t>
    {
        vector<Int_t> layer(nhits);
        for (Int_t i = 0; i < nhits; ++i)
            layer.at(i) = (Int_t) round(Hit_Z.at(i) / Thick);
        return layer;
    }, {"Hit_Z", "nhits"})
    // Number of hits with non-zero energy deposition on each layer
    .Define("hits_on_layer", [] (const vector<Int_t>& layer, const Int_t& nhits)->vector<Int_t>
    {
        vector<Int_t> hits_on_layer(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t ilayer = layer.at(i);
            ++hits_on_layer.at(ilayer);
        }
        return hits_on_layer;
    }, {"layer", "nhits"})
    // Energy deposition on each layer
    .Define("layer_energy", [] (const vector<Int_t>& layer, const vector<Double_t>& Hit_Energy, const Int_t& nhits)->vector<Double_t>
    {
        vector<Double_t> layer_energy(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
            layer_energy.at(layer.at(i)) += Hit_Energy.at(i);
        return layer_energy;
    }, {"layer", "Hit_Energy", "nhits"})
    // Total energy deposition
    .Define("Edep", [] (const vector<Double_t>& Hit_Energy)->Double_t
    {
        Double_t sum = 0;
        for (Double_t i : Hit_Energy)
            sum += i;
        return sum;
    }, {"Hit_Energy"})
    // Average energy deposition of the hits
    .Define("Emean", [] (const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Edep / nhits;
    }, {"Edep", "nhits"})
    // The average centre of gravity, in x direction
    .Define("COG_X_mean", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Energy, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        Double_t cog_x_mean = 0;
        if (Edep == 0)
            return 30.0;
        else
        {
            for (Int_t i = 0; i < nhits; ++i)
                cog_x_mean += Hit_X.at(i) * Hit_Energy.at(i);
            cog_x_mean /= Edep;
            return cog_x_mean;
        }
    }, {"Hit_X", "Hit_Energy", "Edep", "nhits"})
    // The average centre of gravity, in y direction
    .Define("COG_Y_mean", [] (const vector<Double_t>& Hit_Y, const vector<Double_t>& Hit_Energy, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        Double_t cog_y_mean = 0;
        if (Edep == 0)
            return 30.0;
        else
        {
            for (Int_t i = 0; i < nhits; ++i)
                cog_y_mean += Hit_Y.at(i) * Hit_Energy.at(i);
            cog_y_mean /= Edep;
            return cog_y_mean;
        }
    }, {"Hit_Y", "Hit_Energy", "Edep", "nhits"})
    // The average centre of gravity, in z direction
    .Define("COG_Z_mean", [] (const vector<Double_t>& Hit_Z, const vector<Double_t>& Hit_Energy, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        Double_t cog_z_mean = 0;
        if (Edep == 0)
            return 50.0;
        else
        {
            for (Int_t i = 0; i < nhits; ++i)
                cog_z_mean += Hit_Z.at(i) * Hit_Energy.at(i);
            cog_z_mean /= Edep;
            return cog_z_mean;
        }
    }, {"Hit_Z", "Hit_Energy", "Edep", "nhits"})
    // RMS width in x direction (with respect to COGX)
    .Define("xwidth", [] (const vector<Double_t>& Hit_X, const Double_t& COG_X_mean)->Double_t
    {
        TH1D* h1 = new TH1D("h1", "", 120, -30, 30);
        for (Double_t i : Hit_X)
            h1->Fill(i - COG_X_mean);
        const Double_t xwidth = h1->GetRMS();
        delete h1;
        return xwidth;
    }, {"Hit_X", "COG_X_mean"})
    // RMS width in y direction (with respect to COGY)
    .Define("ywidth", [] (const vector<Double_t>& Hit_Y, const Double_t& COG_Y_mean)->Double_t
    {
        TH1D* h1 = new TH1D("h1", "", 120, -30, 30);
        for (Double_t i : Hit_Y)
            h1->Fill(i - COG_Y_mean);
        const Double_t ywidth = h1->GetRMS();
        delete h1;
        return ywidth;
    }, {"Hit_Y", "COG_Y_mean"})
    // RMS depth
    .Define("zdepth", [] (const vector<Double_t>& Hit_Z)->Double_t
    {
        TH1D* h1 = new TH1D("h1", "", nLayer + 1, 0, (nLayer + 1) * Thick);
        for (Double_t i : Hit_Z)
            h1->Fill(i);
        const Double_t zdepth = h1->GetRMS();
        delete h1;
        return zdepth;
    }, {"Hit_Z"})
    // The maximum energy deposition as well as its ID
    .Define("Ecell_max_id", [] (const vector<Int_t>& CellID, const vector<Double_t>& Hit_Energy, const Int_t& nhits)->vector<Double_t>
    {
        vector<Double_t> Ecell_max_id = { 0.0, 0.0 };    // Index 0 for cell ID, index 1 for cell energy
        for (Int_t i = 0; i < nhits; ++i)
            if (Hit_Energy.at(i) > Ecell_max_id.at(1))
            {
                Ecell_max_id.at(0) = CellID.at(i);
                Ecell_max_id.at(1) = Hit_Energy.at(i);
            }
        return Ecell_max_id;
    }, {"CellID", "Hit_Energy", "nhits"})
    .Define("Ecell_max", "Ecell_max_id[1]")
    // The second largest energy deposition as well as its ID
    .Define("Ecell_second_id", [] (const vector<Int_t>& CellID, const vector<Double_t>& Hit_Energy, const Int_t& nhits, const Double_t& Ecell_max)->vector<Double_t>
    {
        vector<Double_t> Ecell_second_id = { 0.0, 0.0 };    // Index 0 for cell ID, index 1 for cell energy
        for (Int_t i = 0; i < nhits; ++i)
            if (Hit_Energy.at(i) > Ecell_second_id.at(1) && Hit_Energy.at(i) < Ecell_max)
            {
                Ecell_second_id.at(0) = CellID.at(i);
                Ecell_second_id.at(1) = Hit_Energy.at(i);
            }
        return Ecell_second_id;
    }, {"CellID", "Hit_Energy", "nhits", "Ecell_max"})
    .Define("Ecell_second", "Ecell_second_id[1]")
    .Define("Emax_sec_diff", "Ecell_max - Ecell_second")
    // Distance between the cell with Ecell_max and the one with Ecell_second
    .Define("Emax_sec_dist", [] (const vector<Double_t>& Ecell_max_id, const vector<Double_t>& Ecell_second_id)->Double_t
    {
        Double_t Emax_sec_dist;
        Double_t x_max = ((Int_t) Ecell_max_id.at(0) % 10000) / 100 * CellWidthX;
        Double_t y_max = (Int_t) Ecell_max_id.at(0) % 100 * CellWidthY;
        const Int_t z_max_temp = (Int_t) Ecell_max_id.at(0) / 10000;
        const Double_t z_max = z_max_temp * Thick;
        if (z_max_temp % 2 == 0)
        {
            x_max += 0.5 * CellWidthX;
            y_max += 0.5 * CellWidthY;
        }
        Double_t x_second = ((Int_t) Ecell_second_id.at(0) % 10000) / 100 * CellWidthX;
        Double_t y_second = (Int_t) Ecell_second_id.at(0) % 100 * CellWidthY;
        const Int_t z_second_temp = (Int_t) Ecell_second_id.at(0) / 10000;
        const Double_t z_second = z_second_temp * Thick;
        if (z_second_temp % 2 == 0)
        {
            x_second += 0.5 * CellWidthX;
            y_second += 0.5 * CellWidthY;
        }
        Emax_sec_dist = Sqrt(Power(x_max - x_second, 2) + Power(y_max - y_second, 2) + Power(z_max - z_second, 2));
        return Emax_sec_dist;
    }, {"Ecell_max_id", "Ecell_second_id"})
    // The energy deposition in the 3*3*3 cells around the one with maximum energy deposition
    .Define("Ecell_max_3", [] (const vector<Int_t>& CellID, const vector<Double_t>& Hit_Energy, const vector<Double_t>& Ecell_max_id)->Double_t
    {
        Double_t Ecell_max_3 = 0.0;
        const Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        const Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        const Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t iz = z - 1; iz <= z + 1; ++iz)
        {
            if (iz < 0 || iz >= nLayer)
                continue;
            for (Int_t ix = x - 1; ix <= x + 1; ++ix)
            {
                if (ix < 0 || ix >= nCellsX)
                    continue;
                for (Int_t iy = y - 1; iy <= y + 1; ++iy)
                {
                    if (iy < 0 || iy >= nCellsY)
                        continue;
                    const Int_t tmp = iz * 10000 + ix * 100 + iy;
                    auto itr = find(CellID.begin(), CellID.end(), tmp);
                    const Int_t index = distance(CellID.begin(), itr);
                    if (index >= CellID.size())
                        continue;
                    Ecell_max_3 += Hit_Energy.at(index);
                }
            }
        }
        return Ecell_max_3;
    }, {"CellID", "Hit_Energy", "Ecell_max_id"})
    // The energy deposition in the 5*5*5 cells around the one with maximum energy deposition
    .Define("Ecell_max_5", [] (const vector<Int_t>& CellID, const vector<Double_t>& Hit_Energy, const vector<Double_t>& Ecell_max_id)->Double_t
    {
        Double_t Ecell_max_5 = 0.0;
        const Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        const Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        const Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t iz = z - 2; iz <= z + 2; ++iz)
        {
            if (iz < 0 || iz >= nLayer)
                continue;
            for (Int_t ix = x - 2; ix <= x + 2; ++ix)
            {
                if (ix < 0 || ix >= nCellsX)
                    continue;
                for (Int_t iy = y - 2; iy <= y + 2; ++iy)
                {
                    if (iy < 0 || iy >= nCellsY)
                        continue;
                    const Int_t tmp = iz * 10000 + ix * 100 + iy;
                    auto itr = find(CellID.begin(), CellID.end(), tmp);
                    const Int_t index = distance(CellID.begin(), itr);
                    if (index >= CellID.size())
                        continue;
                    Ecell_max_5 += Hit_Energy.at(index);
                }
            }
        }
        return Ecell_max_5;
    }, {"CellID", "Hit_Energy", "Ecell_max_id"})
    // The energy deposition in the 7*7*7 cells around the one with maximum energy deposition
    .Define("Ecell_max_7", [] (const vector<Int_t>& CellID, const vector<Double_t>& Hit_Energy, const vector<Double_t>& Ecell_max_id)->Double_t
    {
        Double_t Ecell_max_7 = 0.0;
        const Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        const Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        const Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t iz = z - 3; iz <= z + 3; ++iz)
        {
            if (iz < 0 || iz >= nLayer)
                continue;
            for (Int_t ix = x - 3; ix <= x + 3; ++ix)
            {
                if (ix < 0 || ix >= nCellsX)
                    continue;
                for (Int_t iy = y - 3; iy <= y + 3; ++iy)
                {
                    if (iy < 0 || iy >= nCellsY)
                        continue;
                    const Int_t tmp = iz * 10000 + ix * 100 + iy;
                    auto itr = find(CellID.begin(), CellID.end(), tmp);
                    const Int_t index = distance(CellID.begin(), itr);
                    if (index >= CellID.size())
                        continue;
                    Ecell_max_7 += Hit_Energy.at(index);
                }
            }
        }
        return Ecell_max_7;
    }, {"CellID", "Hit_Energy", "Ecell_max_id"})
    // Energy deposition of the central cell divided by the total energy deposition in the 3*3*3 cells around it
    .Define("E1E3", [] (const Double_t& Ecell_max, const Double_t& Ecell_max_3, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_3;
    }, {"Ecell_max", "Ecell_max_3", "nhits"})
    // Energy deposition of the central cell divided by the total energy deposition in the 5*5*5 cells around it
    .Define("E1E5", [] (const Double_t& Ecell_max, const Double_t& Ecell_max_5, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_5;
    }, {"Ecell_max", "Ecell_max_5", "nhits"})
    // Energy deposition of the central cell divided by the total energy deposition in the 7*7*7 cells around it
    .Define("E1E7", [] (const Double_t& Ecell_max, const Double_t& Ecell_max_7, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_7;
    }, {"Ecell_max", "Ecell_max_7", "nhits"})
    // Energy deposition of the central 3*3*3 cells divided by the total energy deposition in the 5*5*5 cells around it
    .Define("E3E5", [] (const Double_t& Ecell_max_3, const Double_t& Ecell_max_5, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_3 / Ecell_max_5;
    }, {"Ecell_max_3", "Ecell_max_5", "nhits"})
    // Energy deposition of the central 3*3*3 cells divided by the total energy deposition in the 7*7*7 cells around it
    .Define("E3E7", [] (const Double_t& Ecell_max_3, const Double_t& Ecell_max_7, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_3 / Ecell_max_7;
    }, {"Ecell_max_3", "Ecell_max_7", "nhits"})
    // Energy deposition of the cell with largest energy deposition, divided by total energy deposition
    .Define("E1Edep", [] (const Double_t& Ecell_max, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Edep;
    }, {"Ecell_max", "Edep", "nhits"})
    // Energy deposition of the central 3*3*3 cells divided by the total energy deposition
    .Define("E3Edep", [] (const Double_t& Ecell_max_3, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_3 / Edep;
    }, {"Ecell_max_3", "Edep", "nhits"})
    // Energy deposition of the central 5*5*5 cells divided by the total energy deposition
    .Define("E5Edep", [] (const Double_t& Ecell_max_5, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_5 / Edep;
    }, {"Ecell_max_5", "Edep", "nhits"})
    // Energy deposition of the central 7*7*7 cells divided by the total energy deposition
    .Define("E7Edep", [] (const Double_t& Ecell_max_7, const Double_t& Edep, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_7 / Edep;
    }, {"Ecell_max_7", "Edep", "nhits"})
    // RMS value of the positions of all the hits on a layer
    .Define("layer_rms", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Y, const vector<Int_t>& layer, const vector<Double_t>& Hit_Energy, const Int_t& nhits, const Double_t& COG_X_mean, const Double_t& COG_Y_mean)->vector<Double_t>
    {
        vector<Double_t> layer_rms(nLayer);
        vector<TH2D*> hvec;
        hvec.reserve(nLayer);
        for (Int_t i = 0; i < nLayer; ++i)
            hvec.emplace_back(new TH2D("h" + TString(to_string(i)) + "_rms", "Layer RMS", 120, -30, 30, 120, -30, 30));
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t ilayer = layer.at(i);
            hvec.at(ilayer)->Fill(Hit_X.at(i) - COG_X_mean, Hit_Y.at(i) - COG_Y_mean, Hit_Energy.at(i));
        }
        for (Int_t i = 0; i < hvec.size(); ++i)
        {
            if (hvec.at(i)->GetEntries() < 4)
                layer_rms.at(i) = 0.0;
            else
                layer_rms.at(i) = hvec.at(i)->GetRMS();
            delete hvec.at(i);
        }
        vector<TH2D*>().swap(hvec);
        return layer_rms;
    }, {"Hit_X", "Hit_Y", "layer", "Hit_Energy", "nhits", "COG_X_mean", "COG_Y_mean"})
    // The layer where the shower begins; otherwise it is set to be -1
    .Define("shower_start", [] (const vector<Int_t>& hits_on_layer)->Int_t
    {
        Int_t shower_start = -1;
        const Int_t threshold = 3;
        for (Int_t i = 0; i < nLayer - 2; ++i)
            if (hits_on_layer.at(i) >= threshold && hits_on_layer.at(i + 1) >= threshold && hits_on_layer.at(i + 2) >= threshold)
            {
                shower_start = i + 1;
                break;
            }
        return shower_start;
    }, {"hits_on_layer"})
    // The layer where the shower ends
    .Define("shower_end", [] (const vector<Int_t>& hits_on_layer, const Int_t& shower_start)->Int_t
    {
        Int_t shower_end = nLayer;
        const Int_t threshold = 3;
        if (shower_start == -1)
            return -1;
        for (Int_t i = shower_start - 1; i < nLayer - 1; ++i)
            if (hits_on_layer.at(i) < threshold && hits_on_layer.at(i + 1) < threshold)
            {
                shower_end = i + 1;
                break;
            }
        return shower_end;
    }, {"hits_on_layer", "shower_start"})
    // Shower radius with respect to the event axis, which was obtained from 3-dimensional fitting of all the hits
    .Define("shower_radius", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Y, const vector<Double_t>& Hit_Z, const Int_t& nhits)->Double_t
    {
        auto* fit = new Fit3D(Hit_X, Hit_Y, Hit_Z, nhits);
        const Double_t radius = fit->GetRMSRadius();
        delete fit;
        return radius;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    // RMS value of the shower in x direction
    .Define("layer_xwidth", [] (const vector<Double_t>& Hit_X, const vector<Int_t>& layer, const Int_t& nhits, const Double_t& COG_X_mean)->vector<Double_t>
    {
        vector<Double_t> layer_xwidth(nLayer);
        vector<TH1D*> h;
        h.reserve(nLayer);
        for (Int_t l = 0; l < nLayer; ++l)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 120, -30, 30));
        for (Int_t i = 0; i < nhits; ++i)
        {
            const Int_t ilayer = layer.at(i);
            h.at(ilayer)->Fill(Hit_X.at(i) - COG_X_mean);
        }
        for (Int_t i = 0; i < h.size(); ++i)
        {
            layer_xwidth.at(i) = h.at(i)->GetRMS();
            delete h.at(i);
        }
        vector<TH1D*>().swap(h);
        return layer_xwidth;
    }, {"Hit_X", "layer", "nhits", "COG_X_mean"})
    // RMS value of the shower in y direction
    .Define("layer_ywidth", [] (const vector<Double_t>& Hit_Y, const vector<Int_t>& layer, const Int_t& nhits, const Double_t& COG_Y_mean)->vector<Double_t>
    {
        vector<Double_t> layer_ywidth(nLayer);
        vector<TH1D*> h;
        h.reserve(nLayer);
        for (Int_t l = 0; l < nLayer; ++l)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 120, -30, 30));
        for (Int_t i = 0; i < nhits; ++i)
        {
            const Int_t ilayer = layer.at(i);
            h.at(ilayer)->Fill(Hit_Y.at(i) - COG_Y_mean);
        }
        for (Int_t i = 0; i < h.size(); ++i)
        {
            layer_ywidth.at(i) = h.at(i)->GetRMS();
            delete h.at(i);
        }
        vector<TH1D*>().swap(h);
        return layer_ywidth;
    }, {"Hit_Y", "layer", "nhits", "COG_Y_mean"})
    // Number of layers with xwidth, ywidth >= 6 cm
    .Define("shower_layer", [] (const vector<Double_t>& layer_xwidth, const vector<Double_t>& layer_ywidth)->Int_t
    {
        Int_t shower_layer = 0;
        for (Int_t i = 0; i < nLayer; ++i)
            if (layer_xwidth.at(i) >= 6 && layer_ywidth.at(i) >= 6)
                ++shower_layer;
        return shower_layer;
    }, {"layer_xwidth", "layer_ywidth"})
    // Number of layers with at least one hit
    .Define("hit_layer", [] (const vector<Int_t>& layer)->Int_t
    {
        Int_t hit_layer = 0;
        unordered_map<Int_t, Int_t> map_layer_hit;
        for (Int_t i : layer)
            ++map_layer_hit[i];
        for (Int_t i = 0; i < nLayer; ++i)
            if (map_layer_hit.count(i) > 0)
                ++hit_layer;
        return hit_layer;
    }, {"layer"})
    // The proportion of layers with xwidth, ywidth >= 6 cm within the layers with at least one hit
    .Define("shower_layer_ratio", [] (const Int_t& shower_layer, const Int_t& hit_layer, const Int_t& nhits)->Double_t
    {
        if (nhits == 0)
            return 0.0;
        else
            return (Double_t) shower_layer / (Double_t) hit_layer;
    }, {"shower_layer", "hit_layer", "nhits"})
    // Average number of hits in the 3*3*3 cells around a given one
    .Define("shower_density", [] (const vector<Int_t>& CellID, const Int_t& nhits)->Double_t
    {
        Double_t shower_density = 0.0;
        if (nhits == 0)
            return shower_density;
        unordered_map<Int_t, Int_t> map_CellID;
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t index = CellID.at(i);
            ++map_CellID[index];
        }
        for (Int_t j = 0; j < nhits; ++j)
        {
            const Int_t x = (CellID.at(j) % 10000) / 100;
            const Int_t y = CellID.at(j) % 100;
            const Int_t z = CellID.at(j) / 10000;
            for (Int_t iz = z - 1; iz <= z + 1; ++iz)
            {
                if (iz < 0 || iz >= nLayer)
                    continue;
                for (Int_t ix = x - 1; ix <= x + 1; ++ix)
                {
                    if (ix < 0 || ix >= nCellsX)
                        continue;
                    for (Int_t iy = y - 1; iy <= y + 1; ++iy)
                    {
                        if (iy < 0 || iy >= nCellsY)
                            continue;
                        const Int_t tmp = z * 10000 + ix * 100 + iy;
                        shower_density += map_CellID[tmp];
                    }
                }
            }
        }
        shower_density /= nhits;
        return shower_density;
    }, {"CellID", "nhits"})
    // The distance between the layer with largest RMS value of position (with respect to COGX and COGY) and the beginning layer
    .Define("shower_length", [] (const vector<Double_t>& layer_rms, const Int_t& shower_start)->Int_t
    {
        Int_t shower_length;
        Int_t max_layer = 0;
        Double_t max_rms = 0.0;
        for (Int_t i = 0; i < layer_rms.size(); ++i)
            if (layer_rms.at(i) > max_rms)
            {
                max_layer = i + 1;
        	    max_rms = layer_rms.at(i);
            }
//        auto maxPosition = max_element(layer_rms.begin() + shower_start, layer_rms.end());
//        Int_t max_layer = maxPosition - layer_rms.begin();
        if (shower_start >= max_layer || shower_start == -1)
            shower_length = 0;
        else
            shower_length = max_layer - shower_start + 1;
        return shower_length;
    }, {"layer_rms", "shower_start"})
    // 2-dimensional fractal dimension
    .Define("FD_2D", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Y, const vector<Double_t>& Hit_Z, const Int_t& nhits)->vector<Double_t>
    {
        vector<Int_t> scale = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
        const Int_t num = scale.size();
        vector<Double_t> fd_2d(num);
        vector<Int_t> nhits_new(num);
        for (Int_t i = 0; i < num; ++i)
        {
            nhits_new.at(i) = NewScale(Hit_X, Hit_Y, Hit_Z, scale.at(i), scale.at(i), 1);
            if (nhits == 0 || nhits_new.at(i) <= 0)
            {
                fd_2d.at(i) = -1.0;
                continue;
            }
            fd_2d.at(i) = Log((Double_t) nhits / nhits_new.at(i)) / Log(scale.at(i));
        }
        return fd_2d;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    // 3-dimensional fractal dimension
    .Define("FD_3D", [] (const vector<Double_t>& Hit_X, const vector<Double_t>& Hit_Y, const vector<Double_t>& Hit_Z, const Int_t& nhits)->vector<Double_t>
    {
        vector<Int_t> scale = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
        const Int_t num = scale.size();
        vector<Double_t> fd_3d(num);
        vector<Int_t> nhits_new(num);
        for (Int_t i = 0; i < num; ++i)
        {
            nhits_new.at(i) = NewScale(Hit_X, Hit_Y, Hit_Z, scale.at(i), scale.at(i), scale.at(i));
            if (nhits == 0 || nhits_new.at(i) <= 0)
            {
                fd_3d.at(i) = -1.0;
                continue;
            }
            fd_3d.at(i) = Log((Double_t) nhits / nhits_new.at(i)) / Log(scale.at(i));
        }
        return fd_3d;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    // Average value of all the 2-dimensional fractal dimensions
    .Define("FD_2D_mean", [] (const vector<Double_t>& FD_2D)->Double_t
    {
        const Int_t num = FD_2D.size();
        Double_t total = 0;
        for (Int_t i = 0; i < num; ++i)
            total += FD_2D.at(i);
        const Double_t FD_2D_mean = total / num;
        return FD_2D_mean;
    }, {"FD_2D"})
    // Average value of all the 3-dimensional fractal dimensions
    .Define("FD_3D_mean", [] (const vector<Double_t>& FD_3D)->Double_t
    {
        const Int_t num = FD_3D.size();
        Double_t total = 0;
        for (Int_t i = 0; i < num; ++i)
            total += FD_3D.at(i);
        const Double_t FD_3D_mean = total / num;
        return FD_3D_mean;
    }, {"FD_3D"})
    // RMS value of all the 2-dimensional fractal dimensions
    .Define("FD_2D_rms", [] (const vector<Double_t>& FD_2D)->Double_t
    {
        const Int_t num = FD_2D.size();
        Double_t total2 = 0;
        for (Int_t i = 0; i < num; ++i)
            total2 += Power(FD_2D.at(i), 2);
        const Double_t FD_2D_rms = Sqrt(total2 / num);
        return FD_2D_rms;
    }, {"FD_2D"})
    // RMS value of all the 3-dimensional fractal dimensions
    .Define("FD_3D_rms", [] (const vector<Double_t>& FD_3D)->Double_t
    {
        const Int_t num = FD_3D.size();
        Double_t total2 = 0;
        for (Int_t i = 0; i < num; ++i)
            total2 += Power(FD_3D.at(i), 2);
        const Double_t FD_3D_rms = Sqrt(total2 / num);
        return FD_3D_rms;
    }, {"FD_3D"})
    .Snapshot(tree, outname);
    delete dm;

    TFile* f = new TFile((TString) outname, "READ");
    TTree* t = f->Get<TTree>((TString) tree);
    t->SetBranchStatus("*", true);
    const vector<TString> deactivate = { "CellID", "Ecell_max_id", "Ecell_second_id", "FD_2D", "FD_3D", "Hit_Energy", "Hit_X", "Hit_Y", "Hit_Z", "hits_on_layer", "layer", "layer_energy", "layer_rms", "layer_xwidth", "layer_ywidth" };
    for (const TString& de : deactivate)
        t->SetBranchStatus(de, false);
    TFile* fnew = new TFile((TString) outname, "RECREATE");
    TTree* tnew = t->CloneTree();
    tnew->Write(nullptr, TObject::kWriteDelete, 0);
    f->Close();
    fnew->Close();

    return 0;
}
