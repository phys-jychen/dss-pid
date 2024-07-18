#include "PID.h"

Int_t main(Int_t argc, Char_t* argv[])
{
    string file, tree = "dp";
    Int_t sel = 0, hit = 0;
    Int_t rec = 0;
    Int_t bdt = 0, classify = 0;
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

        p->AddVar("COG_X_mean",         'D');
        p->AddVar("COG_Y_mean",         'D');
        p->AddVar("COG_Z_mean",         'D');

        p->AddVar("E1E3",               'D');
        p->AddVar("E1Edep",             'D');
        p->AddVar("E3E5",               'D');
        p->AddVar("E3E7",               'D');
        p->AddVar("E3Edep",             'D');
        p->AddVar("E5Edep",             'D');
        p->AddVar("E7Edep",             'D');
        p->AddVar("ECAL_Cluster_N",     'I');
        p->AddVar("Ecell_max",          'D');
        p->AddVar("Ecell_max_3",        'D');
        p->AddVar("Ecell_max_5",        'D');
        p->AddVar("Ecell_max_7",        'D');
        p->AddVar("Ecell_second",       'D');
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
        p->AddVar("shower_layer_ratio", 'D');
        p->AddVar("shower_length",      'I');
        p->AddVar("shower_radius",      'D');
        p->AddVar("shower_start",       'I');
        p->AddVar("xwidth",             'D');
        p->AddVar("ywidth",             'D');
        p->AddVar("zdepth",             'D');

        const Int_t energy_points = 200;
        string path = "/lustre/collider/chenjiyuan/dss-pid/run/e-signal/root/";

        /* | Particle | Training  |    Test    |
         * | --------------------------------- |
         * |    e-    |  1--200   |  401--600  |
         * |   pi-    | 201--400  |  601--800  |
         * |  gamma   | 801--1000 | 1000--1200 |
         */

        for (Int_t i = 1; i <= energy_points; ++i)
        {
            // Signal
            p->AddTrainSig(path + "training/job" + to_string(i) + "_e-_" + to_string(10 * i) + "MeV/e-_" + to_string(10 * i) + "MeV.root", tree);
            p->AddTestSig( path + "test/job" + to_string(400 + i) + "_e-_" + to_string(10 * i) + "MeV/e-_" + to_string(10 * i) + "MeV.root", tree);

            // Background
            p->AddTrainBkg(path + "training/job" + to_string(200 + i) + "_pi-_" + to_string(10 * i) + "MeV/pi-_" + to_string(10 * i) + "MeV.root", tree);
            p->AddTestBkg( path + "test/job" + to_string(600 + i) + "_pi-_" + to_string(10 * i) + "MeV/pi-_" + to_string(10 * i) + "MeV.root", tree);
//            p->AddTrainBkg(path + "training/job" + to_string(800 + i) + "_gamma_" + to_string(10 * i) + "MeV/gamma_" + to_string(10 * i) + "MeV.root", tree);
//            p->AddTestBkg( path + "test/job" + to_string(1000 + i) + "_gamma_" + to_string(10 * i) + "MeV/gamma_" + to_string(10 * i) + "MeV.root", tree);
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
    else if (help == 0)
    {
        cout << "Invalid input." << endl;
        cout << "Run \"iPID --help\" to display help information." << endl << endl;
    }

    delete p;
    return 0;
}
