#include "PID.h"

Int_t PID::TrainBDT()
{
    TMVA::Tools::Instance();

    unordered_map<TTree*, pair<TString, Double_t>> trsig;
    unordered_map<TTree*, pair<TString, Double_t>> tesig;

    trsig.reserve(train_sig.size());
    tesig.reserve(test_sig.size());

    for (const auto& i : train_sig)
    {
        TFile* f = new TFile(i.first.first, "READ");
        TTree* t = f->Get<TTree>(i.first.second);
        trsig.insert(pair<TTree*, pair<TString, Double_t>>(t, i.second));
    }

    for (const auto& j : test_sig)
    {
        TFile* f = new TFile(j.first.first, "READ");
        TTree* t = f->Get<TTree>(j.first.second);
        tesig.insert(pair<TTree*, pair<TString, Double_t>>(t, j.second));
    }

    // Create a ROOT output file where ntuples, histograms, etc. are stored
    const TString type = "Multiclass";
    const TString outfileName = "TMVA" + type + ".root";
    TFile* outputFile = new TFile( outfileName, "RECREATE" );

    TMVA::Factory* factory = new TMVA::Factory( "TMVA" + type, outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=" + type );
    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    // Variables and spectators should be added here
    for (const auto& i : var)
        dataloader->AddVariable(i.first, i.second);
    for (const auto& j : spec)
        dataloader->AddSpectator(j.first, j.second);

    // Signal and background trees should be added here
    TCut cut = "nhits > 0 && Acts_TagTrk_No == 1 && Acts_RecTrk_No == 1 && Edep_HCAL < 140 && Edep_SideHCAL < 100 && Ecell_max_HCAL < 220 && Ecell_max_SideHCAL < 100";

    for (const auto& i : trsig)
        dataloader->AddTree(i.first, i.second.first, i.second.second, cut, TMVA::Types::kTraining);
    for (const auto& j : tesig)
        dataloader->AddTree(j.first, j.second.first, j.second.second, cut, TMVA::Types::kTesting);

    dataloader->PrepareTrainingAndTestTree( cut, "SplitMode=Random:NormMode=NumEvents:!V" );

    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                         "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2" );

    // Now you can ask the factory to train, test, and evaluate the variables

    // Training and evaluation, using the training and test datasets respectively
    factory->TrainAllMethods();
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // Save the output
    outputFile->Close();

    cout << "==> ROOT file written to " << outputFile->GetName() << endl;
    cout << "==> TMVA" << type << " finished!" << endl;

    delete factory;
    delete dataloader;
    Clear();

    return 0;	
}

Int_t PID::BDTNtuple(const string& file, const string& tree)
{
    EnableImplicitMT();
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "bdt_" + outname;

    TMVA::Tools::Instance();

    const TString type = "Multiclass";
    cout << "==> Start TMVA" << type << "Application" << endl;
    TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );

    Float_t  bdt_Acts_RecTrk_No;
    Float_t  bdt_Acts_TagTrk_No;
    Float_t  bdt_eventid_high;
    Float_t  bdt_eventid_low;
    Float_t  bdt_run_number;

    Float_t  bdt_COG_X_mean;
    Float_t  bdt_COG_Y_mean;
    Float_t  bdt_COG_Z_mean;

    Float_t  bdt_E1E3_centre;
    Float_t  bdt_E1Edep;
    Float_t  bdt_E3E5_centre;
//    Float_t  bdt_E3E7_centre;
    Float_t  bdt_E3Edep;
//    Float_t  bdt_E5Edep;
    Float_t  bdt_E7Edep;
//    Float_t  bdt_ECAL_Cluster_N;
    Float_t  bdt_Ecell_max;
    Float_t  bdt_Ecell_second;
    Float_t  bdt_Ecentre;
    Float_t  bdt_Ecentre_3;
//    Float_t  bdt_Ecentre_5;
    Float_t  bdt_Ecentre_7;
    Float_t  bdt_Eclus_max;
//    Float_t  bdt_Eclus_max_sec_diff;
    Float_t  bdt_Eclus_max_sec_dist;
//    Float_t  bdt_Eclus_second;
    Float_t  bdt_Edep;
//    Float_t  bdt_Emax_sec_diff;
    Float_t  bdt_Emax_sec_dist;
    Float_t  bdt_Emean;

//    Float_t  bdt_FD_2D_mean;
//    Float_t  bdt_FD_2D_rms;
//    Float_t  bdt_FD_3D_mean;
    Float_t  bdt_FD_3D_rms;

    Float_t  bdt_Missing2_Trk_P;

    Float_t  bdt_hit_layer;
    Float_t  bdt_nhits;
    Float_t  bdt_shower_density;
    Float_t  bdt_shower_end;
    Float_t  bdt_shower_layer;
    Float_t  bdt_shower_length;
    Float_t  bdt_shower_radius;
//    Float_t  bdt_weighted_radius;
    Float_t  bdt_xwidth;
    Float_t  bdt_ywidth;
    Float_t  bdt_zdepth;

    reader->AddSpectator("Acts_RecTrk_No", &bdt_Acts_RecTrk_No);
    reader->AddSpectator("Acts_TagTrk_No", &bdt_Acts_TagTrk_No);
    reader->AddSpectator("EventID_High",   &bdt_eventid_high);
    reader->AddSpectator("EventID_Low",    &bdt_eventid_low);
    reader->AddSpectator("RunNumber",      &bdt_run_number);

    reader->AddVariable("COG_X_mean",         &bdt_COG_X_mean);
    reader->AddVariable("COG_Y_mean",         &bdt_COG_Y_mean);
    reader->AddVariable("COG_Z_mean",         &bdt_COG_Z_mean);

    reader->AddVariable("E1E3_centre",        &bdt_E1E3_centre);
    reader->AddVariable("E1Edep",             &bdt_E1Edep);
    reader->AddVariable("E3E5_centre",        &bdt_E3E5_centre);
//    reader->AddVariable("E3E7_centre",        &bdt_E3E7_centre);
    reader->AddVariable("E3Edep",             &bdt_E3Edep);
//    reader->AddVariable("E5Edep",             &bdt_E5Edep);
    reader->AddVariable("E7Edep",             &bdt_E7Edep);
//    reader->AddVariable("ECAL_Cluster_N",     &bdt_ECAL_Cluster_N);
    reader->AddVariable("Ecell_max",          &bdt_Ecell_max);
    reader->AddVariable("Ecell_second",       &bdt_Ecell_second);
    reader->AddVariable("Ecentre",            &bdt_Ecentre);
    reader->AddVariable("Ecentre_3",          &bdt_Ecentre_3);
//    reader->AddVariable("Ecentre_5",          &bdt_Ecentre_5);
    reader->AddVariable("Ecentre_7",          &bdt_Ecentre_7);
    reader->AddVariable("Eclus_max",          &bdt_Eclus_max);
//    reader->AddVariable("Eclus_max_sec_diff", &bdt_Eclus_max_sec_diff);
    reader->AddVariable("Eclus_max_sec_dist", &bdt_Eclus_max_sec_dist);
//    reader->AddVariable("Eclus_second",       &bdt_Eclus_second);
    reader->AddVariable("Edep",               &bdt_Edep);
//    reader->AddVariable("Emax_sec_diff",      &bdt_Emax_sec_diff);
    reader->AddVariable("Emax_sec_dist",      &bdt_Emax_sec_dist);
    reader->AddVariable("Emean",              &bdt_Emean);

//    reader->AddVariable("FD_2D_mean",         &bdt_FD_2D_mean);
//    reader->AddVariable("FD_2D_rms",          &bdt_FD_2D_rms);
//    reader->AddVariable("FD_3D_mean",         &bdt_FD_3D_mean);
    reader->AddVariable("FD_3D_rms",          &bdt_FD_3D_rms);

    reader->AddVariable("Missing2_Trk_P",     &bdt_Missing2_Trk_P);

    reader->AddVariable("hit_layer",          &bdt_hit_layer);
    reader->AddVariable("nhits",              &bdt_nhits);
    reader->AddVariable("shower_density",     &bdt_shower_density);
    reader->AddVariable("shower_end",         &bdt_shower_end);
    reader->AddVariable("shower_layer",       &bdt_shower_layer);
    reader->AddVariable("shower_length",      &bdt_shower_length);
    reader->AddVariable("shower_radius",      &bdt_shower_radius);
//    reader->AddVariable("weighted_radius",    &bdt_weighted_radius);
    reader->AddVariable("xwidth",             &bdt_xwidth);
    reader->AddVariable("ywidth",             &bdt_ywidth);
    reader->AddVariable("zdepth",             &bdt_zdepth);

    reader->BookMVA("BDTG method", "dataset/weights/TMVA" + type + "_BDTG.weights.xml");

    vector<string> rdf_input = {};

    rdf_input.emplace_back("COG_X_mean");
    rdf_input.emplace_back("COG_Y_mean");
    rdf_input.emplace_back("COG_Z_mean");

    rdf_input.emplace_back("E1E3_centre");
    rdf_input.emplace_back("E1Edep");
    rdf_input.emplace_back("E3E5_centre");
//    rdf_input.emplace_back("E3E7_centre");
    rdf_input.emplace_back("E3Edep");
//    rdf_input.emplace_back("E5Edep");
    rdf_input.emplace_back("E7Edep");
//    rdf_input.emplace_back("ECAL_Cluster_N");
    rdf_input.emplace_back("Ecell_max");
    rdf_input.emplace_back("Ecell_second");
    rdf_input.emplace_back("Ecentre");
    rdf_input.emplace_back("Ecentre_3");
//    rdf_input.emplace_back("Ecentre_5");
    rdf_input.emplace_back("Ecentre_7");
    rdf_input.emplace_back("Eclus_max");
//    rdf_input.emplace_back("Eclus_max_sec_diff");
    rdf_input.emplace_back("Eclus_max_sec_dist");
//    rdf_input.emplace_back("Eclus_second");
    rdf_input.emplace_back("Edep");
//    rdf_input.emplace_back("Emax_sec_diff");
    rdf_input.emplace_back("Emax_sec_dist");
    rdf_input.emplace_back("Emean");

//    rdf_input.emplace_back("FD_2D_mean");
//    rdf_input.emplace_back("FD_2D_rms");
//    rdf_input.emplace_back("FD_3D_mean");
    rdf_input.emplace_back("FD_3D_rms");

    rdf_input.emplace_back("Missing2_Trk_P");

    rdf_input.emplace_back("hit_layer");
    rdf_input.emplace_back("nhits");
    rdf_input.emplace_back("shower_density");
    rdf_input.emplace_back("shower_end");
    rdf_input.emplace_back("shower_layer");
    rdf_input.emplace_back("shower_length");
    rdf_input.emplace_back("shower_radius");
//    rdf_input.emplace_back("weighted_radius");
    rdf_input.emplace_back("xwidth");
    rdf_input.emplace_back("ywidth");
    rdf_input.emplace_back("zdepth");

    RDataFrame* df = new RDataFrame(tree, file);
    auto bdtout = df->Define("Likelihood", [&]
        (Double_t COG_X_mean,
         Double_t COG_Y_mean,
         Double_t COG_Z_mean,
         Double_t E1E3_centre,
         Double_t E1Edep,
         Double_t E3E5_centre,
//         Double_t E3E7_centre,
         Double_t E3Edep,
//         Double_t E5Edep,
         Double_t E7Edep,
//         Int_t    ECAL_Cluster_N,
         Double_t Ecell_max,
         Double_t Ecell_second,
         Double_t Ecentre,
         Double_t Ecentre_3,
//         Double_t Ecentre_5,
         Double_t Ecentre_7,
         Double_t Eclus_max,
//         Double_t Eclus_max_sec_diff,
         Double_t Eclus_max_sec_dist,
//         Double_t Eclus_second,
         Double_t Edep,
//         Double_t Emax_sec_diff,
         Double_t Emax_sec_dist,
         Double_t Emean,
//         Double_t FD_2D_mean,
//         Double_t FD_2D_rms,
//         Double_t FD_3D_mean,
         Double_t FD_3D_rms,
         Double_t Missing2_Trk_P,
         Int_t    hit_layer,
         Int_t    nhits,
         Double_t shower_density,
         Int_t    shower_end,
         Int_t    shower_layer,
         Int_t    shower_length,
         Double_t shower_radius,
//         Double_t weighted_radius,
         Double_t xwidth,
         Double_t ywidth,
         Double_t zdepth)->vector<Double_t>
    {
        const Int_t Ntypes = 7;
        vector<Double_t> likelihood(Ntypes);

        bdt_COG_X_mean         = COG_X_mean;
        bdt_COG_Y_mean         = COG_Y_mean;
        bdt_COG_Z_mean         = COG_Z_mean;

        bdt_E1E3_centre        = E1E3_centre;
        bdt_E1Edep             = E1Edep;
        bdt_E3E5_centre        = E3E5_centre;
//        bdt_E3E7_centre        = E3E7_centre;
        bdt_E3Edep             = E3Edep;
//        bdt_E5Edep             = E5Edep;
        bdt_E7Edep             = E7Edep;
//        bdt_ECAL_Cluster_N     = ECAL_Cluster_N;
        bdt_Ecell_max          = Ecell_max;
        bdt_Ecell_second       = Ecell_second;
        bdt_Ecentre            = Ecentre;
        bdt_Ecentre_3          = Ecentre_3;
//        bdt_Ecentre_5          = Ecentre_5;
        bdt_Ecentre_7          = Ecentre_7;
        bdt_Eclus_max          = Eclus_max;
//        bdt_Eclus_max_sec_diff = Eclus_max_sec_diff;
        bdt_Eclus_max_sec_dist = Eclus_max_sec_dist;
//        bdt_Eclus_second       = Eclus_second;
        bdt_Edep               = Edep;
//        bdt_Emax_sec_diff      = Emax_sec_diff;
        bdt_Emax_sec_dist      = Emax_sec_dist;
        bdt_Emean              = Emean;

//        bdt_FD_2D_mean         = FD_2D_mean;
//        bdt_FD_2D_rms          = FD_2D_rms;
//        bdt_FD_3D_mean         = FD_3D_mean;
        bdt_FD_3D_rms          = FD_3D_rms;

        bdt_Missing2_Trk_P     = Missing2_Trk_P;

        bdt_hit_layer          = hit_layer;
        bdt_nhits              = nhits;
        bdt_shower_density     = shower_density;
        bdt_shower_end         = shower_end;
        bdt_shower_layer       = shower_layer;
        bdt_shower_length      = shower_length;
        bdt_shower_radius      = shower_radius;
//        bdt_weighted_radius    = weighted_radius;
        bdt_xwidth             = xwidth;
        bdt_ywidth             = ywidth;
        bdt_zdepth             = zdepth;

        for (Int_t i = 0; i < 7; ++i)
            likelihood.at(i) = reader->EvaluateMulticlass("BDTG method")[i];

        return likelihood;
    }, rdf_input)
    .Define("Signal_Likelihood",    "Likelihood[0]")
    .Define("PNTarget_Likelihood",  "Likelihood[1]")
    .Define("GMMTarget_Likelihood", "Likelihood[2]")
    .Define("GMMECAL_Likelihood",   "Likelihood[3]")
    .Define("ENTarget_Likelihood",  "Likelihood[4]")
    .Define("ENECAL_Likelihood",    "Likelihood[5]")
    .Define("Inclusive_Likelihood", "Likelihood[6]")
    .Snapshot(tree, outname);
    delete df;

    return 0;
}
