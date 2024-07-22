#include "PID.h"

Int_t main(Int_t argc, Char_t* argv[])
{
    string file, tree = "dp";
    Int_t sel = 0, hit = 0;
    Int_t rec = 0;
    Int_t bdt = 0, classify = 0;
    Int_t print = 0;
    Int_t help = 0;

    for (Int_t i = 1; i < argc; i++)
    {
        if (string(argv[i]) == string("--help"))
        {
            help = 1;
            cout << endl;
            cout << "****************************************************************" << endl;
            cout << "                        Help information" << endl;
            cout << "****************************************************************" << endl;
            cout << " Select branches:" << endl;
            cout << "     With default tree \'dp\':     iPID -s -f [file]" << endl;
            cout << "     With other specified tree:  iPID -s -f [file] -t [tree]" << endl;
            cout << " Collect original hit information:" << endl;
            cout << "     With default tree \'dp\':     iPID -h -f [file]" << endl;
            cout << "     With other specified tree:  iPID -h -f [file] -t [tree]" << endl;
            cout << " Reconstruct variables:" << endl;
            cout << "     With default tree \'dp\':     iPID -r -f [file]" << endl;
            cout << "     With other specified tree:  iPID -r -f [file] -t [tree]" << endl;
            cout << " BDT training and testing:" << endl;
            cout << "     With default tree \'dp\':     iPID -b" << endl;
            cout << "     With other specified tree:  iPID -b -t [tree]" << endl;
            cout << " Classification:" << endl;
            cout << "     With default tree \'dp\':     iPID -c -f [file]" << endl;
            cout << "     With other specified tree:  iPID -c -f [file] -t [tree]" << endl;
            cout << " Dumping TP and TN values:       iPID -p [TMVA output file]" << endl;
            cout << "****************************************************************" << endl << endl;
            break;
        }

        else if (string(argv[i]) == string("-f"))
            file = string(argv[i + 1]);

        else if (string(argv[i]) == string("-t"))
            tree = string(argv[i + 1]);

        else if (string(argv[i]) == string("-s"))
            sel = 1;

        else if (string(argv[i]) == string("-h"))
            hit = 1;

        else if (string(argv[i]) == string("-r"))
            rec = 1;

        else if (string(argv[i]) == string("-b"))
            bdt = 1;

        else if (string(argv[i]) == string("-c"))
            classify = 1;

        else if (string(argv[i]) == string("-p"))
            print = 1;
    }

    PID* p = new PID();

    if (sel == 1 && !file.empty())
    {
        cout << "--> Saving branches..." << endl;
        cout << "--> File: " << file << endl;
        cout << "--> Tree: " << tree << endl << endl;

        PID::SaveBranches(file, tree);

        cout << "--> Branch selection finished!" << endl;
    }

    else if (hit == 1 && !file.empty())
    {
        cout << "--> Collecting original hits..." << endl;
        cout << "--> File: " << file << endl;
        cout << "--> Tree: " << tree << endl;

        PID::OriginalHits(file, tree);

        cout << "--> Hit collection finished!" << endl;
    }

    else if (rec == 1 && !file.empty())
    {
        cout << "---> Reconstructing variables..." << endl;
        cout << "---> File: " << file << endl;
        cout << "---> Tree: " << tree << endl;

        PID::GenNtuple(file, tree);

        cout << "---> Variable reconstruction finished!" << endl;
    }

    else if (bdt == 1)
    {
        cout << "----> Training and testing..." << endl;
        cout << "----> Tree: " << tree << endl;

        // Add spectators, which are not used in training, test or evaluation, here
        p->AddSpec("EventNumber",      "Event ID");
        p->AddSpec("RecTrk2_track_No", "Track number in recoil tracker");
        p->AddSpec("RunNumber",        "Run ID");
        p->AddSpec("TagTrk2_track_No", "Track number in tagging tracker");

        // Add variables to be trained, tested and evaluated here
        p->AddVar("COG_X_mean",         'D');
        p->AddVar("COG_Y_mean",         'D');
        p->AddVar("COG_Z_mean",         'D');

        p->AddVar("E1E3_centre",        'D');
        p->AddVar("E1Edep",             'D');
        p->AddVar("E3E5_centre",        'D');
        p->AddVar("E3E7_centre",        'D');
        p->AddVar("E3Edep",             'D');
        p->AddVar("E5Edep",             'D');
        p->AddVar("E7Edep",             'D');
        p->AddVar("ECAL_Cluster_N",     'I');
        p->AddVar("Ecell_max",          'D');
        p->AddVar("Ecell_second",       'D');
        p->AddVar("Ecentre",            'D');
        p->AddVar("Ecentre_3",          'D');
        p->AddVar("Ecentre_5",          'D');
        p->AddVar("Ecentre_7",          'D');
        p->AddVar("Eclus_max",          'D');
        p->AddVar("Eclus_max_sec_diff", 'D');
        p->AddVar("Eclus_max_sec_dist", 'D');
        p->AddVar("Eclus_second",       'D');
        p->AddVar("Edep",               'D');
        p->AddVar("Emax_sec_diff",      'D');
        p->AddVar("Emax_sec_dist",      'D');
        p->AddVar("Emean",              'D');

        p->AddVar("FD_2D_mean",         'D');
        p->AddVar("FD_2D_rms",          'D');
        p->AddVar("FD_3D_mean",         'D');
        p->AddVar("FD_3D_rms",          'D');

        p->AddVar("hit_layer",          'I');
        p->AddVar("nhits",              'I');
        p->AddVar("shower_density",     'D');
        p->AddVar("shower_end",         'I');
        p->AddVar("shower_layer",       'I');
        p->AddVar("shower_length",      'I');
        p->AddVar("shower_radius",      'D');
        p->AddVar("weighted_radius",    'D');
        p->AddVar("xwidth",             'D');
        p->AddVar("ywidth",             'D');
        p->AddVar("zdepth",             'D');

        // Add training and test events here
        // Signals: (1, 5, 10, 50, 100, 500, 800, 1000) MeV
        // Backgrounds:  en_ecal, en_target, gmm_ecal, gmm_target, pn_target;  inclusive
        const unordered_map<string, Int_t> bkg_num = { {"en_ecal", 200}, {"en_target", 20}, {"gmm_ecal", 25}, {"gmm_target", 250}, {"pn_target", 15}, {"inclusive", 5} };
        const string path = "/lustre/collider/chenjiyuan/dss-pid/run/dp-signal/";
        const string bkg = "en_ecal";
        const string bkg_path = (bkg == "inclusive") ? "inclusive/root/" : "rare/" + bkg + "/";

        const Int_t mass = 10;    // In MeV
        const Int_t bkg_points = bkg_num.at(bkg);

        p->AddTrainSig(path + "signal/root/training/Mass" + to_string(mass) + "MeV/Mass" + to_string(mass) + "MeV.root", tree);
        p->AddTestSig( path + "signal/root/test/Mass"     + to_string(mass) + "MeV/Mass" + to_string(mass) + "MeV.root", tree);

        for (Int_t i = 1; i <= bkg_points; ++i)
        {
            p->AddTrainBkg(path + bkg_path + "training/job" + to_string(i) + "/" + bkg + "_" + to_string(i) + ".root",                       tree);
            p->AddTestBkg( path + bkg_path + "test/job" + to_string(i + bkg_points) + "/" + bkg + "_" + to_string(i + bkg_points) + ".root", tree);
        }

        p->TrainBDT();

        cout << "----> Training and testing finished!" << endl;
	}

    else if (classify == 1 && !file.empty())
    {
        cout << "----> Classifying..." << endl;
        cout << "----> File: " << file << endl;
        cout << "----> Tree: " << tree << endl << endl;

        PID::BDTNtuple(file, tree);

        cout << "----> Classification finished!" << endl;
    }

    else if (print == 1 && !file.empty())
    {
        cout << "-----> Printing to CSV..." << endl;
        cout << "-----> File: " << file << endl;

        PID::PrintCSV(file);

        cout << "-----> Printing to CSV finished!" << endl;
    }

    else if (help == 0)
    {
        cout << "Invalid input." << endl;
        cout << "Run \"iPID --help\" to display help information." << endl << endl;
    }

    delete p;
    return 0;
}
