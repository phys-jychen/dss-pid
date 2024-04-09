#include "PID.h"

Int_t PID::TrainBDT()
{
    TMVA::Tools::Instance();

    vector<TTree*> trsig;
    vector<TTree*> trbkg;
    vector<TTree*> tesig;
    vector<TTree*> tebkg;

    trsig.reserve(train_sig.size());
    trbkg.reserve(train_bkg.size());
    tesig.reserve(test_sig.size());
    tebkg.reserve(test_bkg.size());

    for (const auto& i : train_sig)
    {
        TFile* f = TFile::Open(i.first, "READ");
        TTree* t = f->Get<TTree>(i.second);
        trsig.emplace_back(t);
    }

    for (const auto& j : train_bkg)
    {
        TFile* f = TFile::Open(j.first, "READ");
        TTree* t = f->Get<TTree>(j.second);
        trbkg.emplace_back(t);
    }

    for (const auto& k : test_sig)
    {
        TFile* f = TFile::Open(k.first, "READ");
        TTree* t = f->Get<TTree>(k.second);
        tesig.emplace_back(t);
    }

    for (const auto& l : test_bkg)
    {
        TFile* f = TFile::Open(l.first, "READ");
        TTree* t = f->Get<TTree>(l.second);
        tebkg.emplace_back(t);
    }

    // Create a ROOT output file where ntuples, histograms, etc. are stored
    TString outfileName( "TMVAMulticlass.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory* factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
    for (const auto& i : var)
        dataloader->AddVariable(i.first, i.second);

    // Signal and background trees should be added here
    for (auto i : trsig)
        dataloader->AddSignalTree(i, 1.0, TMVA::Types::kTraining);
    for (auto j : trbkg)
        dataloader->AddBackgroundTree(j, 1.0, TMVA::Types::kTraining);
    for (auto k : tesig)
        dataloader->AddSignalTree(k, 1.0, TMVA::Types::kTesting);
    for (auto l : tebkg)
        dataloader->AddBackgroundTree(l, 1.0, TMVA::Types::kTesting);

    // Add different cuts for signal and background events here
    TCut cut_sig = "nhits > 0";
    TCut cut_bkg = "nhits > 0";

    dataloader->PrepareTrainingAndTestTree( cut_sig, cut_bkg, "SplitMode=Random:NormMode=NumEvents:!V" );

    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2" );

    // Now you can ask the factory to train, test, and evaluate the variables

    // Training, using the set of training events
    factory->TrainAllMethods();

    // Evaluation, using the set of test events
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    outputFile->Close();

    cout << "==> ROOT file written: " << outputFile->GetName() << endl;
    cout << "==> TMVAClassification finished!" << endl;

    delete factory;
    delete dataloader;
    Clear();

    return 0;	
}

Int_t PID::BDTNtuple(const string& fname, const string& tname)
{
    EnableImplicitMT();
    string outname = fname;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "bdt_" + outname;

    // Loading the library
    TMVA::Tools::Instance();

    // Default MVA methods to be trained + tested
    map<string, Int_t> Use;

    // Cut optimisation
    Use["BDTG"] = 1;
    cout << "==> Start TMVAMulticlassApplication" << endl;
    TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );

    Float_t  bdt_COG_X_mean;
    Float_t  bdt_COG_Y_mean;
    Float_t  bdt_COG_Z_mean;

    Float_t  bdt_E1E9;
    Float_t  bdt_E1Edep;
    Float_t  bdt_E25Edep;
    Float_t  bdt_E49Edep;
    Float_t  bdt_E9E25;
    Float_t  bdt_E9E49;
    Float_t  bdt_E9Edep;
    Float_t  bdt_ECAL_Cluster_N;
    Float_t  bdt_Ecell_max;
    Float_t  bdt_Ecell_max_25;
    Float_t  bdt_Ecell_max_49;
    Float_t  bdt_Ecell_max_9;
    Float_t  bdt_Ecell_second;
    Float_t  bdt_Eclus_max;
    Float_t  bdt_Eclus_max_sec_diff;
    Float_t  bdt_Eclus_max_sec_dist;
    Float_t  bdt_Eclus_second;
    Float_t  bdt_Edep;
    Float_t  bdt_Emax_sec_diff;
    Float_t  bdt_Emax_sec_dist;
    Float_t  bdt_Emean;

    Float_t  bdt_FD_2D_mean;
    Float_t  bdt_FD_2D_rms;
    Float_t  bdt_FD_3D_mean;
    Float_t  bdt_FD_3D_rms;

    Float_t  bdt_hit_layer;
    Float_t  bdt_nhits;
    Float_t  bdt_shower_density;
    Float_t  bdt_shower_end;
    Float_t  bdt_shower_layer;
    Float_t  bdt_shower_layer_ratio;
    Float_t  bdt_shower_length;
    Float_t  bdt_shower_radius;
    Float_t  bdt_shower_start;
    Float_t  bdt_xwidth;
    Float_t  bdt_ywidth;
    Float_t  bdt_zdepth;

    reader->AddVariable("COG_X_mean",         &bdt_COG_X_mean);
    reader->AddVariable("COG_Y_mean",         &bdt_COG_Y_mean);
    reader->AddVariable("COG_Z_mean",         &bdt_COG_Z_mean);

    reader->AddVariable("E1E9",               &bdt_E1E9);
    reader->AddVariable("E1Edep",             &bdt_E1Edep);
    reader->AddVariable("E25Edep",            &bdt_E25Edep);
    reader->AddVariable("E49Edep",            &bdt_E49Edep);
    reader->AddVariable("E9E25",              &bdt_E9E25);
    reader->AddVariable("E9E49",              &bdt_E9E49);
    reader->AddVariable("E9Edep",             &bdt_E9Edep);
    reader->AddVariable("ECAL_Cluster_N",     &bdt_ECAL_Cluster_N);
    reader->AddVariable("Ecell_max",          &bdt_Ecell_max);
    reader->AddVariable("Ecell_max_25",       &bdt_Ecell_max_25);
    reader->AddVariable("Ecell_max_49",       &bdt_Ecell_max_49);
    reader->AddVariable("Ecell_max_9",        &bdt_Ecell_max_9);
    reader->AddVariable("Ecell_second",       &bdt_Ecell_second);
    reader->AddVariable("Eclus_max",          &bdt_Eclus_max);
    reader->AddVariable("Eclus_max_sec_diff", &bdt_Eclus_max_sec_diff);
    reader->AddVariable("Eclus_max_sec_dist", &bdt_Eclus_max_sec_dist);
    reader->AddVariable("Eclus_second",       &bdt_Eclus_second);
    reader->AddVariable("Edep",               &bdt_Edep);
    reader->AddVariable("Emax_sec_diff",      &bdt_Emax_sec_diff);
    reader->AddVariable("Emax_sec_dist",      &bdt_Emax_sec_dist);
    reader->AddVariable("Emean",              &bdt_Emean);

    reader->AddVariable("FD_2D_mean",         &bdt_FD_2D_mean);
    reader->AddVariable("FD_2D_rms",          &bdt_FD_2D_rms);
    reader->AddVariable("FD_3D_mean",         &bdt_FD_3D_mean);
    reader->AddVariable("FD_3D_rms",          &bdt_FD_3D_rms);

    reader->AddVariable("hit_layer",          &bdt_hit_layer);
    reader->AddVariable("nhits",              &bdt_nhits);
    reader->AddVariable("shower_density",     &bdt_shower_density);
    reader->AddVariable("shower_end",         &bdt_shower_end);
    reader->AddVariable("shower_layer",       &bdt_shower_layer);
    reader->AddVariable("shower_layer_ratio", &bdt_shower_layer_ratio);
    reader->AddVariable("shower_length",      &bdt_shower_length);
    reader->AddVariable("shower_radius",      &bdt_shower_radius);
    reader->AddVariable("shower_start",       &bdt_shower_start);
    reader->AddVariable("xwidth",             &bdt_xwidth);
    reader->AddVariable("ywidth",             &bdt_ywidth);
    reader->AddVariable("zdepth",             &bdt_zdepth);

    reader->BookMVA("BDTG method", TString("dataset/weights/TMVAMulticlass_BDTG.weights.xml"));
    cout << "Booked" << endl;

    vector<string> rdf_input = {};

    rdf_input.emplace_back("COG_X_mean");
    rdf_input.emplace_back("COG_Y_mean");
    rdf_input.emplace_back("COG_Z_mean");

    rdf_input.emplace_back("E1E9");
    rdf_input.emplace_back("E1Edep");
    rdf_input.emplace_back("E25Edep");
    rdf_input.emplace_back("E49Edep");
    rdf_input.emplace_back("E9E25");
    rdf_input.emplace_back("E9E49");
    rdf_input.emplace_back("E9Edep");
    rdf_input.emplace_back("ECAL_Cluster_N");
    rdf_input.emplace_back("Ecell_max");
    rdf_input.emplace_back("Ecell_max_25");
    rdf_input.emplace_back("Ecell_max_49");
    rdf_input.emplace_back("Ecell_max_9");
    rdf_input.emplace_back("Ecell_second");
    rdf_input.emplace_back("Eclus_max");
    rdf_input.emplace_back("Eclus_max_sec_diff");
    rdf_input.emplace_back("Eclus_max_sec_dist");
    rdf_input.emplace_back("Eclus_second");
    rdf_input.emplace_back("Edep");
    rdf_input.emplace_back("Emax_sec_diff");
    rdf_input.emplace_back("Emax_sec_dist");
    rdf_input.emplace_back("Emean");

    rdf_input.emplace_back("FD_2D_mean");
    rdf_input.emplace_back("FD_2D_rms");
    rdf_input.emplace_back("FD_3D_mean");
    rdf_input.emplace_back("FD_3D_rms");

    rdf_input.emplace_back("hit_layer");
    rdf_input.emplace_back("nhits");
    rdf_input.emplace_back("shower_density");
    rdf_input.emplace_back("shower_end");
    rdf_input.emplace_back("shower_layer");
    rdf_input.emplace_back("shower_layer_ratio");
    rdf_input.emplace_back("shower_length");
    rdf_input.emplace_back("shower_radius");
    rdf_input.emplace_back("shower_start");
    rdf_input.emplace_back("xwidth");
    rdf_input.emplace_back("ywidth");
    rdf_input.emplace_back("zdepth");

    RDataFrame df(tname, fname);

    auto bdtout = df.Define("Response_Sig", [&]
        (Double_t COG_X_mean,
         Double_t COG_Y_mean,
         Double_t COG_Z_mean,
         Double_t E1E9,
         Double_t E1Edep,
         Double_t E25Edep,
         Double_t E49Edep,
         Double_t E9E25,
         Double_t E9E49,
         Double_t E9Edep,
         Int_t    ECAL_Cluster_N,
         Double_t Ecell_max,
         Double_t Ecell_max_25,
         Double_t Ecell_max_49,
         Double_t Ecell_max_9,
         Double_t Ecell_second,
         Double_t Eclus_max,
         Double_t Eclus_max_sec_diff,
         Double_t Eclus_max_sec_dist,
         Double_t Eclus_second,
         Double_t Edep,
         Double_t Emax_sec_diff,
         Double_t Emax_sec_dist,
         Double_t Emean,
         Double_t FD_2D_mean,
         Double_t FD_2D_rms,
         Double_t FD_3D_mean,
         Double_t FD_3D_rms,
         Int_t    hit_layer,
         Int_t    nhits,
         Double_t shower_density,
         Int_t    shower_end,
         Int_t    shower_layer,
         Double_t shower_layer_ratio,
         Int_t    shower_length,
         Double_t shower_radius,
         Int_t    shower_start,
         Double_t xwidth,
         Double_t ywidth,
         Double_t zdepth)->Float_t
    {
        bdt_COG_X_mean         = COG_X_mean;
        bdt_COG_Y_mean         = COG_Y_mean;
        bdt_COG_Z_mean         = COG_Z_mean;

        bdt_E1E9               = E1E9;
        bdt_E1Edep             = E1Edep;
        bdt_E25Edep            = E25Edep;
        bdt_E49Edep            = E49Edep;
        bdt_E9E25              = E9E25;
        bdt_E9E49              = E9E49;
        bdt_E9Edep             = E9Edep;
        bdt_ECAL_Cluster_N     = ECAL_Cluster_N;
        bdt_Ecell_max          = Ecell_max;
        bdt_Ecell_max_25       = Ecell_max_25;
        bdt_Ecell_max_49       = Ecell_max_49;
        bdt_Ecell_max_9        = Ecell_max_9;
        bdt_Ecell_second       = Ecell_second;
        bdt_Eclus_max          = Eclus_max;
        bdt_Eclus_max_sec_diff = Eclus_max_sec_diff;
        bdt_Eclus_max_sec_dist = Eclus_max_sec_dist;
        bdt_Eclus_second       = Eclus_second;
        bdt_Edep               = Edep;
        bdt_Emax_sec_diff      = Emax_sec_diff;
        bdt_Emax_sec_dist      = Emax_sec_dist;
        bdt_Emean              = Emean;

        bdt_FD_2D_mean         = FD_2D_mean;
        bdt_FD_2D_rms          = FD_2D_rms;
        bdt_FD_3D_mean         = FD_3D_mean;
        bdt_FD_3D_rms          = FD_3D_rms;

        bdt_hit_layer          = hit_layer;
        bdt_nhits              = nhits;
        bdt_shower_density     = shower_density;
        bdt_shower_end         = shower_end;
        bdt_shower_layer       = shower_layer;
        bdt_shower_layer_ratio = shower_layer_ratio;
        bdt_shower_length      = shower_length;
        bdt_shower_radius      = shower_radius;
        bdt_shower_start       = shower_start;
        bdt_xwidth             = xwidth;
        bdt_ywidth             = ywidth;
        bdt_zdepth             = zdepth;
        return (reader->EvaluateMulticlass( "BDTG method" ))[0];
    }, rdf_input)
    .Snapshot(tname, outname);
    return 0;
}
