#include "PID.h"

Int_t PID::PrintCSV(const string& file)
{
    TFile* f = new TFile((TString) file, "READ");
    TDirectoryFile* d = f->Get<TDirectoryFile>("dataset");
    d->cd();
    TTree* t = d->Get<TTree>("TestTree");

    const Long64_t signal_total = t->GetEntries("classID == 0");
    const Long64_t pn_target_total = t->GetEntries("classID == 1");
    const Long64_t gmm_target_total = t->GetEntries("classID == 2");
    const Long64_t gmm_ecal_total = t->GetEntries("classID == 3");
    const Long64_t en_target_total = t->GetEntries("classID == 4");
    const Long64_t en_ecal_total = t->GetEntries("classID == 5");
    const Long64_t inclusive_total = t->GetEntries("classID == 6");

    const Int_t nsteps = 1000000;

    ofstream num("./Numbers.csv");
    num << "threshold,"
        << "signal_total,nsig_signal,nbkg_signal,"
        << "pn_target_total,nsig_pn_target,nbkg_pn_target,"
        << "gmm_target_total,nsig_gmm_target,nbkg_gmm_target,"
        << "gmm_ecal_total,nsig_gmm_ecal,nbkg_gmm_ecal,"
        << "en_target_total,nsig_en_target,nbkg_en_target,"
        << "en_ecal_total,nsig_en_ecal,nbkg_en_ecal,"
        << "inclusive_total,nsig_inclusive,nbkg_inclusive\n";

    cout << "Cycle begins!" << endl;

    for (Int_t i = 0; i < nsteps; ++i)
    {
        if (i > 0.001 * nsteps && i < 0.999 * nsteps && (i + 1) % 1000 != 0)
            continue;

        const Double_t selection = (Double_t) i / nsteps;

        Long64_t nsig_signal = t->GetEntries("classID == 0 && signal > " + (TString) to_string(selection));
        Long64_t nbkg_signal = t->GetEntries("classID != 0 && signal > " + (TString) to_string(selection));
        Long64_t nsig_pn_target = t->GetEntries("classID == 1 && pn_target > " + (TString) to_string(selection));
        Long64_t nbkg_pn_target = t->GetEntries("classID != 1 && pn_target > " + (TString) to_string(selection));
        Long64_t nsig_gmm_target = t->GetEntries("classID == 2 && gmm_target > " + (TString) to_string(selection));
        Long64_t nbkg_gmm_target = t->GetEntries("classID != 2 && gmm_target > " + (TString) to_string(selection));
        Long64_t nsig_gmm_ecal = t->GetEntries("classID == 3 && gmm_ecal > " + (TString) to_string(selection));
        Long64_t nbkg_gmm_ecal = t->GetEntries("classID != 3 && gmm_ecal > " + (TString) to_string(selection));
        Long64_t nsig_en_target = t->GetEntries("classID == 4 && en_target > " + (TString) to_string(selection));
        Long64_t nbkg_en_target = t->GetEntries("classID != 4 && en_target > " + (TString) to_string(selection));
        Long64_t nsig_en_ecal = t->GetEntries("classID == 5 && en_ecal > " + (TString) to_string(selection));
        Long64_t nbkg_en_ecal = t->GetEntries("classID != 5 && en_ecal > " + (TString) to_string(selection));
        Long64_t nsig_inclusive = t->GetEntries("classID == 6 && inclusive > " + (TString) to_string(selection));
        Long64_t nbkg_inclusive = t->GetEntries("classID != 6 && inclusive > " + (TString) to_string(selection));

        num << selection << ","
            << signal_total << "," << nsig_signal << "," << nbkg_signal << ","
            << pn_target_total << "," << nsig_pn_target << "," << nbkg_pn_target << ","
            << gmm_target_total << "," << nsig_gmm_target << "," << nbkg_gmm_target << ","
            << gmm_ecal_total << "," << nsig_gmm_ecal << "," << nbkg_gmm_ecal << ","
            << en_target_total << "," << nsig_en_target << "," << nbkg_en_target << ","
            << en_ecal_total << "," << nsig_en_ecal << "," << nbkg_en_ecal << ","
            << inclusive_total << "," << nsig_inclusive << "," << nbkg_inclusive << "\n";

        if (((i <= 0.001 * nsteps || i >= 0.999 * nsteps) && (i + 1) % 100 == 0)
            || (i > 0.001 * nsteps && i < 0.999 * nsteps && (i + 1) % 20000 == 0))
            cout << "Cycle " << i + 1 << " / " << nsteps << " finished!" << endl;
    }

    num.close();
    f->Close();

    return 0;
}