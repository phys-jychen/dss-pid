#include "Variables.h"

const Double_t CellWidthX = 2.5, CellWidthY = 2.5;
const Double_t Thick = 4.0;

const Int_t nCellX = 21, nCellY = 21;
const Int_t nLayer = 11;

const Double_t BiasX = 0.5 * (nCellX - 1) * CellWidthX;
const Double_t BiasY = 0.5 * (nCellY - 1) * CellWidthY;

Int_t NHScaleV2(vector<Double_t> pos_x, vector<Double_t> pos_y, vector<Double_t> pos_z, Int_t RatioX, Int_t RatioY, Int_t RatioZ)
{
    Int_t ReScaledNH = 0;
    Int_t tmpI = 0;
    Int_t tmpJ = 0;
    Int_t tmpK = 0;
    Double_t tmpEn = 0;
    Int_t NewCellID0 = 0;
    const Int_t NumHit = pos_x.size();

    map<Double_t, Double_t> testIDtoEnergy;

    for (Int_t i = 0; i < NumHit; ++i)
    {
        Double_t x = pos_x.at(i);
        Double_t y = pos_y.at(i);
        Double_t z = pos_z.at(i);

        tmpI = ((Int_t) ((x + BiasX) / CellWidthX) + (Int_t) (Abs(x) / x)) / RatioX;
        tmpJ = ((Int_t) ((y + BiasY) / CellWidthY) + (Int_t) (Abs(y) / y)) / RatioY;
        tmpK = (Int_t) (z / Thick) / RatioZ;
        tmpEn = 1;

        NewCellID0 = (tmpK << 24) + (tmpJ << 12) + tmpI;

        if (testIDtoEnergy.find(NewCellID0) == testIDtoEnergy.end())
            testIDtoEnergy[NewCellID0] = tmpEn;
        else
            testIDtoEnergy[NewCellID0] += tmpEn;
    }

    ReScaledNH = testIDtoEnergy.size();
    return ReScaledNH;
}

Variables::Variables() {}

Variables::~Variables() {}

Int_t Variables::GenNtuple(const string& file, const string& tree)
{
//    EnableImplicitMT();
    DisableImplicitMT();
    RDataFrame* dm = new RDataFrame(tree, file);
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "rec_" + outname;
    auto fout = dm->Define("nhits", "(Int_t) Hit_X.size()")
    // Transfer z position to the layer number
    .Define("layer", [] (vector<Double_t> Hit_Z, Int_t nhits)
    {
        vector<Int_t> layer = {};
        for (Int_t i = 0; i < nhits; ++i)
            layer.emplace_back((Int_t) round(Hit_Z.at(i) / Thick));
        return layer;
    }, {"Hit_Z", "nhits"})
    // Number of hits with non-zero energy deposition on each layer
    .Define("hits_on_layer", [] (vector<Int_t> layer, vector<Double_t> Hit_Energy, Int_t nhits)
    {
        vector<Int_t> hits_on_layer(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t ilayer = layer.at(i);
            ++hits_on_layer.at(ilayer);
        }
        return hits_on_layer;
    }, {"layer", "Hit_Energy", "nhits"})
    // Energy deposition on each layer
    .Define("layer_energy", [] (vector<Int_t> layer, vector<Double_t> Hit_Energy, Int_t nhits)
    {
        vector<Double_t> layer_energy(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
            layer_energy.at(layer.at(i)) += Hit_Energy.at(i);
        return layer_energy;
    }, {"layer", "Hit_Energy", "nhits"})
    // Total energy deposition
    .Define("Edep", [] (vector<Double_t> Hit_Energy)
    {
        Double_t sum = 0;
        for (Double_t i : Hit_Energy)
            sum += i;
        return sum;
    }, {"Hit_Energy"})
    // Average energy deposition of the hits
    .Define("Emean", [] (Double_t Edep, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Edep / nhits;
    }, {"Edep", "nhits"})
    /*
    // Centre of gravity of each layer, in x direction
    .Define("COG_X", [] (vector<Double_t> Hit_X, vector<Int_t> layer, vector<Double_t> Hit_Energy, vector<Double_t> layer_energy, Int_t nhits)
    {
        vector<Double_t> cog_x(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
            cog_x.at(layer.at(i)) += Hit_X.at(i) * Hit_Energy.at(i);
        for (Int_t j = 0; j < cog_x.size(); ++j)
        {
            if (cog_x.at(j) == 0)
                continue;
            else
                cog_x.at(j) /= layer_energy.at(j);
        }
        return cog_x;
    }, {"Hit_X", "layer", "Hit_Energy", "layer_energy", "nhits"})
    // Centre of gravity of each layer, in y direction
    .Define("COG_Y", [] (vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Hit_Energy, vector<Double_t> layer_energy, Int_t nhits)
    {
        vector<Double_t> cog_y(nLayer);
        for (Int_t i = 0; i < nhits; ++i)
            cog_y.at(layer.at(i)) += Hit_Y.at(i) * Hit_Energy.at(i);
        for (Int_t j = 0; j < cog_y.size(); ++j)
        {
            if (cog_y.at(j) == 0)
                continue;
            else
                cog_y.at(j) /= layer_energy.at(j);
        }
        return cog_y;
    }, {"Hit_Y", "layer", "Hit_Energy", "layer_energy", "nhits"})
    // Centre of gravity of every 4 layers, in x direction
    .Define("COG_X_4", [] (vector<Double_t> Hit_X, vector<Int_t> layer, vector<Double_t> Hit_Energy, vector<Double_t> layer_energy, Int_t nhits)
    {
        vector<Double_t> cog_x_4(nLayer / 4);
        vector<Double_t> energy_4layer(nLayer / 4);
        for (Int_t i = 0; i < nhits; ++i)
            cog_x_5.at(layer.at(i) / 4) += Hit_X.at(i) * Hit_Energy.at(i);
        for (Int_t j = 0; j < nLayer; ++j)
            energy_5layer.at(j / 4) += layer_energy.at(j);
        for (Int_t k = 0; k < nLayer / 4; ++k)
        {
            if (cog_x_4.at(k) == 0)
                cog_x_4.at(k) = 30.0;
            else
                cog_x_4.at(k) /= energy_4layer.at(k);
        }
        return cog_x_4;
    }, {"Hit_X", "layer", "Hit_Energy", "layer_energy", "nhits"})
    .Define("COG_X_4_0", "COG_X_4[0]")
    .Define("COG_X_4_1", "COG_X_4[1]")
    .Define("COG_X_4_2", "COG_X_4[2]")
    // Centre of gravity of every 5 layers, in y direction
    .Define("COG_Y_4", [] (vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Hit_Energy, vector<Double_t> layer_energy, Int_t nhits)
    {
        vector<Double_t> cog_y_4(nLayer / 4);
        vector<Double_t> energy_4layer(nLayer / 4);
        for (Int_t i = 0; i < nhits; ++i)
            cog_y_5.at(layer.at(i) / 4) += Hit_Y.at(i) * Hit_Energy.at(i);
        for (Int_t j = 0; j < nLayer; ++j)
            energy_5layer.at(j / 4) += layer_energy.at(j);
        for (Int_t k = 0; k < nLayer / 4; ++k)
        {
            if (cog_y_4.at(k) == 0)
                cog_y_4.at(k) = 30.0;
            else
                cog_y_4.at(k) /= energy_4layer.at(k);
        }
        return cog_y_4;
    }, {"Hit_Y", "layer", "Hit_Energy", "layer_energy", "nhits"})
    .Define("COG_Y_4_0", "COG_Y_4[0]")
    .Define("COG_Y_4_1", "COG_Y_4[1]")
    .Define("COG_Y_4_2", "COG_Y_4[2]")
    // Centre of gravity of every 5 layers, in z direction
    .Define("COG_Z_4", [] (vector<Double_t> Hit_Z, vector<Int_t> layer, vector<Double_t> Hit_Energy, vector<Double_t> layer_energy, Int_t nhits)
    {
        vector<Double_t> cog_z_4(nLayer / 4);
        vector<Double_t> energy_4layer(nLayer / 4);
        for (Int_t i = 0; i < nhits; ++i)
            cog_z_5.at(layer.at(i) / 4) += Hit_Z.at(i) * Hit_Energy.at(i);
        for (Int_t j = 0; j < nLayer; ++j)
            energy_5layer.at(j / 4) += layer_energy.at(j);
        for (Int_t k = 0; k < nLayer / 4; ++k)
        {
            if (cog_z_4.at(k) == 0)
                cog_z_4.at(k) = 50.0;
            else
                cog_z_4.at(k) /= energy_4layer.at(k);
        }
        return cog_z_4;
    }, {"Hit_Z", "layer", "Hit_Energy", "layer_energy", "nhits"})
    .Define("COG_Z_4_0", "COG_Z_4[0]")
    .Define("COG_Z_4_1", "COG_Z_4[1]")
    .Define("COG_Z_4_2", "COG_Z_4[2]")
    */
    // The average centre of gravity, in x direction
    .Define("COG_X_mean", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Energy, Double_t Edep, Int_t nhits)
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
    .Define("COG_Y_mean", [] (vector<Double_t> Hit_Y, vector<Double_t> Hit_Energy, Double_t Edep, Int_t nhits)
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
    .Define("COG_Z_mean", [] (vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy, Double_t Edep, Int_t nhits)
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
    .Define("xwidth", [] (vector<Double_t> Hit_X, Double_t COG_X_mean)
    {
        TH1D* h1 = new TH1D("h1", "", 100, -30, 30);
        for (Double_t i : Hit_X)
            h1->Fill(i - COG_X_mean);
        Double_t xwidth = h1->GetRMS();
        delete h1;
        return xwidth;
    }, {"Hit_X", "COG_X_mean"})
    // RMS width in y direction (with respect to COGY)
    .Define("ywidth", [] (vector<Double_t> Hit_Y, Double_t COG_Y_mean)
    {
        TH1D* h1 = new TH1D("h1", "", 100, -30, 30);
        for (Double_t i : Hit_Y)
            h1->Fill(i - COG_Y_mean);
        Double_t ywidth = h1->GetRMS();
        delete h1;
        return ywidth;
    }, {"Hit_Y", "COG_Y_mean"})
    // RMS depth
    .Define("zdepth", [] (vector<Double_t> Hit_Z)
    {
        TH1D* h1 = new TH1D("h1", "", nLayer + 1, 0, (nLayer + 1) * Thick);
        for (Double_t i : Hit_Z)
            h1->Fill(i);
        Double_t zdepth = h1->GetRMS();
        delete h1;
        return zdepth;
    }, {"Hit_Z"})
    // The energy deposition of the fired cells
    .Define("Ecell", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Hit_Energy, Int_t nhits)
    {
        unordered_map<Int_t, Double_t> Ecell_map;
        vector<vector<Double_t>> Ecell = { {}, {} };
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t x = round((Hit_X.at(i) + BiasX) / CellWidthX);
            Int_t y = round((Hit_Y.at(i) + BiasY) / CellWidthY);
            Int_t index = 10000 * layer.at(i) + 100 * x + y;
            Ecell_map[index] += Hit_Energy.at(i);
        }
        for (auto it = Ecell_map.cbegin(); it != Ecell_map.cend(); ++it)
        {
            Ecell.at(0).emplace_back(it->first);
            Ecell.at(1).emplace_back(it->second);
        }
        return Ecell;
    }, {"Hit_X", "Hit_Y", "layer", "Hit_Energy", "nhits"})
    // The maximum energy deposition as well as the ID of the fired cells
    .Define("Ecell_max_id", [] (vector<vector<Double_t>> Ecell, Int_t nhits)
    {
        vector<Double_t> Ecell_max_id = { (Double_t) nhits, 0.0 };
        for (Int_t i = 0; i < nhits; ++i)
            if (Ecell.at(1).at(i) > Ecell_max_id.at(1))
            {
                Ecell_max_id.at(0) = Ecell.at(0).at(i);
                Ecell_max_id.at(1) = Ecell.at(1).at(i);
            }
        return Ecell_max_id;
    }, {"Ecell", "nhits"})
    .Define("Ecell_max", "Ecell_max_id[1]")
    // The energy deposition in the 3*3 cells around the one with maximum energy deposition
    .Define("Ecell_max_9", [] (vector<vector<Double_t>> Ecell, vector<Double_t> Ecell_max_id)
    {
        Double_t Ecell_max_9 = 0.0;
        Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t ix = x - 1; ix <= x + 1; ++ix)
        {
            if (ix < 0 || ix >= nCellX)
                continue;
            for (Int_t iy = y - 1; iy <= y + 1; ++iy)
            {
                if (iy < 0 || iy >= nCellY)
                    continue;
                Int_t tmp = z * 10000 + ix * 100 + iy;
                auto itr = find(Ecell.at(0).begin(), Ecell.at(0).end(), tmp);
                Int_t index = distance(Ecell.at(0).begin(), itr);
                if (index >= Ecell.at(0).size())
                    continue;
                Ecell_max_9 += Ecell.at(1).at(index);
            }
        }
        return Ecell_max_9;
    }, {"Ecell", "Ecell_max_id"})
    // The energy deposition in the 5*5 cells around the one with maximum energy deposition
    .Define("Ecell_max_25", [] (vector<vector<Double_t>> Ecell, vector<Double_t> Ecell_max_id)
    {
        Double_t Ecell_max_25 = 0.0;
        Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t ix = x - 2; ix <= x + 2; ++ix)
        {
            if (ix < 0 || ix >= nCellX)
                continue;
            for (Int_t iy = y - 2; iy <= y + 2; ++iy)
            {
                if (iy < 0 || iy >= nCellY)
                    continue;
                Int_t tmp = z * 10000 + ix * 100 + iy;
                auto itr = find(Ecell.at(0).begin(), Ecell.at(0).end(), tmp);
                Int_t index = distance(Ecell.at(0).begin(), itr);
                if (index >= Ecell.at(0).size())
                    continue;
                Ecell_max_25 += Ecell.at(1).at(index);
            }
        }
        return Ecell_max_25;
    }, {"Ecell", "Ecell_max_id"})
    // The energy deposition in the 7*7 cells around the one with maximum energy deposition
    .Define("Ecell_max_49", [] (vector<vector<Double_t>> Ecell, vector<Double_t> Ecell_max_id)
    {
        Double_t Ecell_max_49 = 0.0;
        Int_t z = (Int_t) Ecell_max_id.at(0) / 10000;
        Int_t x = ((Int_t) Ecell_max_id.at(0) % 10000) / 100;
        Int_t y = (Int_t) Ecell_max_id.at(0) % 100;
        for (Int_t ix = x - 3; ix <= x + 3; ++ix)
        {
            if (ix < 0 || ix >= nCellX)
                continue;
            for (Int_t iy = y - 3; iy <= y + 3; ++iy)
            {
                if (iy < 0 || iy >= nCellY)
                    continue;
                Int_t tmp = z * 10000 + ix * 100 + iy;
                auto itr = find(Ecell.at(0).begin(), Ecell.at(0).end(), tmp);
                Int_t index = distance(Ecell.at(0).begin(), itr);
                if (index >= Ecell.at(0).size())
                    continue;
                Ecell_max_49 += Ecell.at(1).at(index);
            }
        }
        return Ecell_max_49;
    }, {"Ecell", "Ecell_max_id"})
    // Energy deposition of the central cell divided by the total energy deposition in the 3*3 cells around it
    .Define("E1E9", [] (Double_t Ecell_max, Double_t Ecell_max_9, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_9;
    }, {"Ecell_max", "Ecell_max_9", "nhits"})
    // Energy deposition of the central cell divided by the total energy deposition in the 5*5 cells around it
    .Define("E1E25", [] (Double_t Ecell_max, Double_t Ecell_max_25, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_25;
    }, {"Ecell_max", "Ecell_max_25", "nhits"})
    // Energy deposition of the central cell divided by the total energy deposition in the 7*7 cells around it
    .Define("E1E49", [] (Double_t Ecell_max, Double_t Ecell_max_49, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max / Ecell_max_49;
    }, {"Ecell_max", "Ecell_max_49", "nhits"})
    // Energy deposition of the central 3*3 cells divided by the total energy deposition in the 5*5 cells around it
    .Define("E9E25", [] (Double_t Ecell_max_9, Double_t Ecell_max_25, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_9 / Ecell_max_25;
    }, {"Ecell_max_9", "Ecell_max_25", "nhits"})
    // Energy deposition of the central 3*3 cells divided by the total energy deposition in the 7*7 cells around it
    .Define("E9E49", [] (Double_t Ecell_max_9, Double_t Ecell_max_49, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return Ecell_max_9 / Ecell_max_49;
    }, {"Ecell_max_9", "Ecell_max_49", "nhits"})
    // RMS value of the positions of all the hits on a layer
    .Define("layer_rms", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Hit_Energy, Int_t nhits, Double_t COG_X_mean, Double_t COG_Y_mean)->vector<Double_t>
    {
        vector<Double_t> layer_rms(nLayer);
        vector<TH2D*> hvec;
        for (Int_t i = 0; i < nLayer; ++i)
            hvec.emplace_back(new TH2D("h" + TString(to_string(i)) + "_rms", "Layer RMS", 100, -30, 30, 100, -30, 30));
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
    // The layer where the shower begins; otherwise it is set to be nLayer + 2
    .Define("shower_start", [] (vector<Int_t> hits_on_layer)
    {
        Int_t shower_start = nLayer + 2;
        const Int_t threshold = 4;
        for (Int_t i = 0; i < nLayer - 3; ++i)
            if (hits_on_layer.at(i) >= threshold && hits_on_layer.at(i + 1) >= threshold && hits_on_layer.at(i + 2) >= threshold && hits_on_layer.at(i + 3) >= threshold)
            {
                shower_start = i;
                break;
            }
        return shower_start;
    }, {"hits_on_layer"})
    // The layer where the shower ends
    .Define("shower_end", [] (vector<Int_t> hits_on_layer, Int_t shower_start)
    {
        Int_t shower_end = nLayer + 2;
        const Int_t threshold = 4;
        if (shower_start == nLayer + 2)
            return shower_end;
        for (Int_t i = shower_start; i < hits_on_layer.size() - 2; ++i)
            if (hits_on_layer.at(i) < threshold && hits_on_layer.at(i + 1) < threshold && hits_on_layer.at(i + 2) < threshold)
            {
                shower_end = i;
                break;
            }
        if (shower_end == nLayer + 2)
            shower_end = nLayer;
        return shower_end;
    }, {"hits_on_layer", "shower_start"})
    // Shower radius with respect to COGX and COGY (between beginning and ending layers)
    .Define("shower_radius", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, Int_t beginning, Int_t ending, Int_t nhits, Double_t COG_X_mean, Double_t COG_Y_mean)
    {
        Double_t d2 = 0;
        Int_t hits = 0;
        for (Int_t i = 0; i < nhits; ++i)
            if (layer.at(i) >= beginning && layer.at(i) < ending)
            {
                ++hits;
                d2 += Power(Hit_X.at(i) - COG_X_mean, 2) + Power(Hit_Y.at(i) - COG_Y_mean, 2);
            }
        if (hits == 0)
            return 0.0;
        Double_t radius = Sqrt(d2 / hits);
        return radius;
    }, {"Hit_X", "Hit_Y", "layer", "shower_start", "shower_end", "nhits", "COG_X_mean", "COG_Y_mean"})
    // RMS value of the shower in x direction
    .Define("layer_xwidth", [] (vector<Double_t> Hit_X, vector<Int_t> layer, Int_t nhits, Double_t COG_X_mean)
    {
        vector<Double_t> layer_xwidth(nLayer);
        vector<TH1D*> h;
        for (Int_t l = 0; l < nLayer; ++l)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -30, 30));
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t ilayer = layer.at(i);
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
    .Define("layer_ywidth", [] (vector<Double_t> Hit_Y, vector<Int_t> layer, Int_t nhits, Double_t COG_Y_mean)
    {
        vector<Double_t> layer_ywidth(nLayer);
        vector<TH1D*> h;
        for (Int_t l = 0; l < nLayer; ++l)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -30, 30));
        for (Int_t i = 0; i < nhits; ++i)
        {
            Int_t ilayer = layer.at(i);
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
    // Number of layers with xwidth, ywidth >= 60 mm
    .Define("shower_layer", [] (vector<Double_t> layer_xwidth, vector<Double_t> layer_ywidth)
    {
        Int_t shower_layer = 0;
        for (Int_t i = 0; i < nLayer; ++i)
            if (layer_xwidth.at(i) >= 6 && layer_ywidth.at(i) >= 6)
                ++shower_layer;
        return shower_layer;
    }, {"layer_xwidth", "layer_ywidth"})
    // Number of layers with at least one hit
    .Define("hit_layer", [] (vector<Int_t> layer)
    {
        Int_t hit_layer = 0;
        unordered_map<Int_t, Int_t> map_layer_hit;
        for (Double_t i : layer)
            ++map_layer_hit[i];
        for (Int_t i = 0; i < nLayer; ++i)
            if (map_layer_hit.count(i) > 0)
                ++hit_layer;
        return hit_layer;
    }, {"layer"})
    // The proportion of layers with xwidth, ywidth >= 60 mm within the layers with at least one hit
    .Define("shower_layer_ratio", [] (Int_t shower_layer, Int_t hit_layer, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        else
            return (Double_t) shower_layer / (Double_t) hit_layer;
    }, {"shower_layer", "hit_layer", "nhits"})
    // Average number of hits in the 3*3 cells around a given one
    .Define("shower_density", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Hit_Energy, Int_t nhits)
    {
        Double_t shower_density = 0.0;
        if (nhits == 0)
            return shower_density;
        unordered_map<Int_t, Int_t> map_CellID;
        for (Int_t j = 0; j < nhits; ++j)
        {
            Int_t x = round((Hit_X.at(j) + BiasX) / CellWidthX);
            Int_t y = round((Hit_Y.at(j) + BiasY) / CellWidthY);
            Int_t z = layer.at(j);
            Int_t index = z * 10000 + x * 100 + y;
            map_CellID[index] += 1;
        }
        for (Int_t i = 0; i < nhits; ++i)
        {
            if (Hit_Energy.at(i) == 0.0)
                continue;
            Int_t x = round((Hit_X.at(i) + BiasX) / CellWidthX);
            Int_t y = round((Hit_Y.at(i) + BiasY) / CellWidthY);
            Int_t z = layer.at(i);
            for (Int_t ix = x - 1; ix <= x + 1; ++ix)
            {
                if (ix < 0 || ix > nCellX - 1)
                    continue;
                for (Int_t iy = y - 1; iy <= y + 1; ++iy)
                {
                    if (iy < 0 || iy > nCellY - 1)
                        continue;
                    Int_t tmp = z * 10000 + ix * 100 + iy;
                    shower_density += map_CellID[tmp];
                }
            }
        }
        shower_density /= nhits;
        return shower_density;
    }, {"Hit_X", "Hit_Y", "layer", "Hit_Energy", "nhits"})
    // The distance between the layer with largest RMS value of position (with respect to COGX and COGY) and the beginning layer
    .Define("shower_length", [] (vector<Double_t> layer_rms, Int_t shower_start)
    {
        Int_t shower_length = 0;
        Int_t max_layer = 0;
        Double_t max_rms = 0.0;
        for (Int_t i = 0; i < layer_rms.size(); ++i)
            if (layer_rms.at(i) > max_rms)
            {
                max_layer = i;
        	    max_rms = layer_rms.at(i);
            }
//        auto maxPosition = max_element(layer_rms.begin() + shower_start, layer_rms.end());
//        Int_t max_layer = maxPosition - layer_rms.begin();
        if (shower_start >= max_layer)
            shower_length = 0;
        else
            shower_length = max_layer - shower_start;
        return shower_length;
    }, {"layer_rms", "shower_start"})
    // 2-dimensional fractal dimension
    .Define("FD_2D", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, Int_t nhits)
    {
        vector<Int_t> scale = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150 };
        const Int_t num = scale.size();
        vector<Double_t> fd_2d(num);
        vector<Int_t> NResizeHit(num);
        for (Int_t i = 0; i < num; ++i)
        {
            NResizeHit.at(i) = NHScaleV2(Hit_X, Hit_Y, Hit_Z, scale.at(i), scale.at(i), 1);
            if (nhits == 0 || NResizeHit.at(i) <= 0)
            {
                fd_2d.at(i) = -1.0;
                continue;
            }
            fd_2d.at(i) = Log((Double_t) nhits / NResizeHit.at(i)) / Log((Double_t) scale.at(i));
        }
        return fd_2d;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    .Define("FD_2D_2",   "FD_2D[0]")
    .Define("FD_2D_3",   "FD_2D[1]")
    .Define("FD_2D_4",   "FD_2D[2]")
    .Define("FD_2D_5",   "FD_2D[3]")
    .Define("FD_2D_6",   "FD_2D[4]")
    .Define("FD_2D_7",   "FD_2D[5]")
    .Define("FD_2D_8",   "FD_2D[6]")
    .Define("FD_2D_9",   "FD_2D[7]")
    .Define("FD_2D_10",  "FD_2D[8]")
    .Define("FD_2D_12",  "FD_2D[9]")
    .Define("FD_2D_15",  "FD_2D[10]")
    .Define("FD_2D_20",  "FD_2D[11]")
    .Define("FD_2D_30",  "FD_2D[12]")
    .Define("FD_2D_40",  "FD_2D[13]")
    .Define("FD_2D_50",  "FD_2D[14]")
    .Define("FD_2D_60",  "FD_2D[15]")
    .Define("FD_2D_70",  "FD_2D[16]")
    .Define("FD_2D_80",  "FD_2D[17]")
    .Define("FD_2D_90",  "FD_2D[18]")
    .Define("FD_2D_100", "FD_2D[19]")
    .Define("FD_2D_110", "FD_2D[20]")
    .Define("FD_2D_120", "FD_2D[21]")
    .Define("FD_2D_130", "FD_2D[22]")
    .Define("FD_2D_140", "FD_2D[23]")
    .Define("FD_2D_150", "FD_2D[24]")
    // 3-dimensional fractal dimension
    .Define("FD_3D", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, Int_t nhits)
    {
        vector<Int_t> scale = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150 };
        const Int_t num = scale.size();
        vector<Double_t> fd_3d(num);
        vector<Int_t> NResizeHit(num);
        for (Int_t i = 0; i < num; ++i)
        {
            NResizeHit.at(i) = NHScaleV2(Hit_X, Hit_Y, Hit_Z, scale.at(i), scale.at(i), scale.at(i));
            if (nhits == 0 || NResizeHit.at(i) <= 0)
            {
                fd_3d.at(i) = -1.0;
                continue;
            }
            fd_3d.at(i) = Log((Double_t) nhits / NResizeHit.at(i)) / Log((Double_t) scale.at(i));
        }
        return fd_3d;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    .Define("FD_3D_2",   "FD_3D[0]")
    .Define("FD_3D_3",   "FD_3D[1]")
    .Define("FD_3D_4",   "FD_3D[2]")
    .Define("FD_3D_5",   "FD_3D[3]")
    .Define("FD_3D_6",   "FD_3D[4]")
    .Define("FD_3D_7",   "FD_3D[5]")
    .Define("FD_3D_8",   "FD_3D[6]")
    .Define("FD_3D_9",   "FD_3D[7]")
    .Define("FD_3D_10",  "FD_3D[8]")
    .Define("FD_3D_12",  "FD_3D[9]")
    .Define("FD_3D_15",  "FD_3D[10]")
    .Define("FD_3D_20",  "FD_3D[11]")
    .Define("FD_3D_30",  "FD_3D[12]")
    .Define("FD_3D_40",  "FD_3D[13]")
    .Define("FD_3D_50",  "FD_3D[14]")
    .Define("FD_3D_60",  "FD_3D[15]")
    .Define("FD_3D_70",  "FD_3D[16]")
    .Define("FD_3D_80",  "FD_3D[17]")
    .Define("FD_3D_90",  "FD_3D[18]")
    .Define("FD_3D_100", "FD_3D[19]")
    .Define("FD_3D_110", "FD_3D[20]")
    .Define("FD_3D_120", "FD_3D[21]")
    .Define("FD_3D_130", "FD_3D[22]")
    .Define("FD_3D_140", "FD_3D[23]")
    .Define("FD_3D_150", "FD_3D[24]")
    // Average value of all the 2-dimensional fractal dimensions
    .Define("FD_2D_mean", [] (vector<Double_t> FD_2D)
    {
        const Int_t num = FD_2D.size();
        Double_t total = 0;
        for (Int_t i = 0; i < num; ++i)
            total += FD_2D.at(i);
        Double_t FD_2D_mean = total / num;
        return FD_2D_mean;
    }, {"FD_2D"})
    // Average value of all the 3-dimensional fractal dimensions
    .Define("FD_3D_mean", [] (vector<Double_t> FD_3D)
    {
        const Int_t num = FD_3D.size();
        Double_t total = 0;
        for (Int_t i = 0; i < num; ++i)
            total += FD_3D.at(i);
        Double_t FD_3D_mean = total / num;
        return FD_3D_mean;
    }, {"FD_3D"})
    // RMS value of all the 2-dimensional fractal dimensions
    .Define("FD_2D_rms", [] (vector<Double_t> FD_2D)
    {
        const Int_t num = FD_2D.size();
        Double_t total2 = 0;
        for (Int_t i = 0; i < num; ++i)
            total2 += Power(FD_2D.at(i), 2);
        Double_t FD_2D_rms = Sqrt(total2 / num);
        return FD_2D_rms;
    }, {"FD_2D"})
    // RMS value of all the 3-dimensional fractal dimensions
    .Define("FD_3D_rms", [] (vector<Double_t> FD_3D)
    {
        const Int_t num = FD_3D.size();
        Double_t total2 = 0;
        for (Int_t i = 0; i < num; ++i)
            total2 += Power(FD_3D.at(i), 2);
        Double_t FD_3D_rms = Sqrt(total2 / num);
        return FD_3D_rms;
    }, {"FD_3D"})
    // The number of tracks of an event, with Hough transformation applied
    .Define("ntrack", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
    {
        Int_t ntrack = 0;
        Hough* hough = new Hough(Hit_X, Hit_Y, Hit_Z, Hit_Energy);
        ntrack = hough->GetNtrack();
        delete hough;
        return ntrack;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    /*
    //
    .Define("hclx", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
    {
        vector<Double_t> hclx;
        Hough* hough = new Hough(Hit_X, Hit_Y, Hit_Z, Hit_Energy);
        hclx = hough->GetHclX();
        delete hough;
        return hclx;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    // 
    .Define("hcly", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
    {
        vector<Double_t> hcly;
        Hough* hough = new Hough(Hit_X, Hit_Y, Hit_Z, Hit_Energy);
        hcly = hough->GetHclY();
        delete hough;
        return hcly;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    // 
    .Define("hclz", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
    {
        vector<Double_t> hclz;
        Hough* hough = new Hough(Hit_X, Hit_Y, Hit_Z, Hit_Energy);
        hclz = hough->GetHclZ();
        delete hough;
        return hclz;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    //
    .Define("hcle", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
    {
        vector<Double_t> hcle;
        Hough* hough = new Hough(Hit_X, Hit_Y, Hit_Z, Hit_Energy);
        hcle = hough->GetHclE();
        delete hough;
        return hcle;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    */
//    .Range(1)
    .Snapshot(tree, outname);
    delete dm;

    TFile* f = new TFile((TString) outname, "READ");
    TTree* t = f->Get<TTree>((TString) tree);
    t->SetBranchStatus("*", 1);
    vector<TString> deactivate = { "Ecell", "Ecell_max_id", "FD_2D", "FD_3D", "Hit_Energy", "Hit_Phi", "Hit_Theta", "Hit_X", "Hit_Y", "Hit_Z", "hits_on_layer", "layer", "layer_energy", "layer_rms", "layer_xwidth", "layer_ywidth" };
    for (TString de : deactivate)
        t->SetBranchStatus(de, 0);
    TFile* fnew = new TFile((TString) outname, "RECREATE");
    TTree* tnew = t->CloneTree();
    tnew->Write(0, TObject::kWriteDelete, 0);
    f->Close();
    fnew->Close();

    return 0;
}
