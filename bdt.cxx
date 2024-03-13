#include "BDT.h"
using namespace std;

Int_t main(Int_t argc, Char_t* argv[])
{
    string file = "", tree = "dp";
    Int_t train = 0, bdt = 0, help = 0;

    for (Int_t i = 1; i < argc; i++)
    {
        if (string(argv[i]) == string("-h") || string(argv[i]) == string("-help"))
        {
            help = 1;
            cout << endl;
            cout << "Help information" << endl << endl;
            cout << "Select branches:" << endl;
            cout << "    With default tree \"dp\":     iHit -s -f [file]" << endl;
            cout << "    With other specified tree:  iHit -s -f [file] -t [tree]" << endl;
            cout << "Collect original hit information:" << endl;
            cout << "    With default tree \"dp\":     iHit -f [file]" << endl;
            cout << "    With other specified tree:  iHit -f [file] -t [tree]" << endl;
            cout << "Reconstruct variables:" << endl;
            cout << "    With default tree \"dp\":     iRec -f [file]" << endl;
            cout << "    With other specified tree:  iRec -f [file] -t [tree]" << endl;
            cout << "PID with BDT:" << endl;
            cout << "    With default tree \"dp\":     iBDT -r" << endl;
            cout << "    With other specified tree:  iBDT -r -t [tree]" << endl;
            cout << "Classification:" << endl;
            cout << "    With default tree \"dp\":     iBDT -v -f [file]" << endl;
            cout << "    With other specified tree:  iBDT -v -f [file] -t [tree]" << endl << endl;
            break;
        }

        else if (string(argv[i]) == string("-f"))
            file = string(argv[i + 1]);

        else if (string(argv[i]) == string("-t"))
            tree = string(argv[i + 1]);

        else if (string(argv[i]) == string("-r"))
            train = 1;

        else if (string(argv[i]) == string("-v"))
            bdt = 1;
    }

    BDT* b = new BDT();

    if (train == 1)
    {
        cout << "----> Training and testing..." << endl;
        cout << "----> Tree: " << tree << endl;

        b->AddVar("COG_X_mean",         'D');
        b->AddVar("COG_Y_mean",         'D');
        b->AddVar("COG_Z_mean",         'D');

        b->AddVar("E1E9",               'D');
//        b->AddVar("E1Edep",             'D');
//        b->AddVar("E25Edep",            'D');
//        b->AddVar("E49Edep",            'D');
//        b->AddVar("E9E25",              'D');
//        b->AddVar("E9E49",              'D');
//        b->AddVar("E9Edep",             'D');
//        b->AddVar("ECAL_Cluster_N",     'I');
//        b->AddVar("Ecell_max",          'D');
//        b->AddVar("Ecell_max_25",       'D');
//        b->AddVar("Ecell_max_49",       'D');
//        b->AddVar("Ecell_max_9",        'D');
        b->AddVar("Ecell_second",       'D');
        b->AddVar("Eclus_max",          'D');
//        b->AddVar("Eclus_max_sec_diff", 'D');
//        b->AddVar("Eclus_max_sec_dist", 'D');
//        b->AddVar("Eclus_second",       'D');
        b->AddVar("Edep",               'D');
//        b->AddVar("Emax_sec_diff",      'D');
//        b->AddVar("Emax_sec_dist",      'D');
        b->AddVar("Emean",              'D');

        b->AddVar("FD_2D_mean",         'D');
//        b->AddVar("FD_2D_rms",          'D');
//        b->AddVar("FD_3D_mean",         'D');
//        b->AddVar("FD_3D_rms",          'D');

//        b->AddVar("hit_layer",          'I');
//        b->AddVar("nhits",              'I');
//        b->AddVar("ntrack",             'I');
        b->AddVar("shower_density",     'D');
//        b->AddVar("shower_end",         'I');
//        b->AddVar("shower_layer",       'I');
//        b->AddVar("shower_layer_ratio", 'D');
//        b->AddVar("shower_length",      'I');
//        b->AddVar("shower_radius",      'D');
//        b->AddVar("shower_start",       'I');
//        b->AddVar("xwidth",             'D');
//        b->AddVar("ywidth",             'D');
//        b->AddVar("zdepth",             'D');

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
            b->AddTrainSig(path + "training/job" + to_string(i) + "_e-_" + to_string(10 * i) + "MeV/e-_" + to_string(10 * i) + "MeV.root", tree);
            b->AddTestSig( path + "test/job" + to_string(400 + i) + "_e-_" + to_string(10 * i) + "MeV/e-_" + to_string(10 * i) + "MeV.root", tree);

            // Background
//            b->AddTrainBkg(path + "training/job" + to_string(800 + i) + "_gamma_" + to_string(10 * i) + "MeV/gamma_" + to_string(10 * i) + "MeV.root", tree);
//            b->AddTestBkg( path + "test/job" + to_string(1000 + i) + "_gamma_" + to_string(10 * i) + "MeV/gamma_" + to_string(10 * i) + "MeV.root", tree);
            b->AddTrainBkg(path + "training/job" + to_string(200 + i) + "_pi-_" + to_string(10 * i) + "MeV/pi-_" + to_string(10 * i) + "MeV.root", tree);
            b->AddTestBkg( path + "test/job" + to_string(600 + i) + "_pi-_" + to_string(10 * i) + "MeV/pi-_" + to_string(10 * i) + "MeV.root", tree);
        }

        b->TrainBDT();

        cout << "----> Training and testing finished!" << endl;
	}

    else if (bdt == 1 && file != "")
    {
        cout << "----> Classifying..." << endl;
        cout << "----> File: " << file << endl;
        cout << "----> Tree: " << tree << endl << endl;

        b->BDTNtuple(file, tree);

        cout << "----> Classification finished!" << endl;
    }

    else if (help == 0)
    {
        cout << "Invalid input." << endl;
        cout << "Run \"iBDT -h[elp]\" to display help information." << endl << endl;
    }

    delete b;
    return 0;
}
