#include "PID.h"

Int_t PID::PrintCSV(const string& file)
{
    TFile* f = new TFile((TString) file, "READ");
    TDirectoryFile* d = f->Get<TDirectoryFile>("dataset");
    d->cd();
    TTree* t = d->Get<TTree>("TestTree");

    const Int_t sig_total = t->GetEntries("classID == 0");
    const Int_t bkg_total = t->GetEntries("classID == 1");
    const Int_t nsteps = 1000000;

    ofstream num("./Numbers.csv");
    num << "threshold,"
        << "sig_total,"
        << "bkg_total,"
        << "nsig,"
        << "nbkg\n";

    cout << "Cycle begins!" << endl;

    for (Int_t i = 0; i < nsteps; ++i)
    {
        if (i > 0.001 * nsteps && i < 0.999 * nsteps && (i + 1) % 1000 != 0)
            continue;

        Double_t selection = (Double_t) i / nsteps;
        Int_t nsig = t->GetEntries("classID == 0 && BDTG > " + (TString) to_string(selection));
        Int_t nbkg = t->GetEntries("classID == 1 && BDTG > " + (TString) to_string(selection));

        num << selection << ","
            << sig_total << ","
            << bkg_total << ","
            << nsig << ","
            << nbkg << "\n";

        if (((i <= 0.001 * nsteps || i >= 0.999 * nsteps) && (i + 1) % 100 == 0)
            || (i > 0.001 * nsteps && i < 0.999 * nsteps && (i + 1) % 20000 == 0))
            cout << "Cycle " << i + 1 << " / " << nsteps << " finished!" << endl;
    }

    num.close();
    f->Close();

    return 0;
}