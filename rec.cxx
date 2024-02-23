#include "Variables.h"
using namespace std;

Int_t main(Int_t argc, Char_t* argv[])
{
    string file = "", tree = "dp";
    Int_t help = 0;

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
    }

    Variables* v = new Variables();

    if (file != "")
    {
        cout << "---> Reconstructing variables..." << endl;
        cout << "---> File: " << file << endl;
        cout << "---> Tree: " << tree << endl << endl;

        v->GenNtuple(file, tree);

        cout << "---> Variable reconstruction finished!" << endl;
    }

    else if (help == 0)
    {
        cout << "Invalid input." << endl;
        cout << "Run \"iRec -h[elp]\" to display help information." << endl << endl;
    }

    delete v;
    return 0;
}
