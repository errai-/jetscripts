#include "Test.h"


inline bool compatibility(double mass_sum, double mass_diff) {
    return (mass_sum < 400 && mass_sum > 300 && mass_diff < 50);
}

inline bool mass_study(double m1, double m2, double n1, double n2, bool noisy) {
    double sum_1 = m1+n2, sum_2 = m2+n1;
    double diff_1 = fabs(m1-n2), diff_2 = fabs(m2-n1);
    
    unsigned success_count = 0;
    if (compatibility(sum_1,diff_1)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m1 << "Jet t " << n2 << endl;
    }
    if (compatibility(sum_2,diff_2)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m2 << " Jet t " << n1 << endl;
    }
    
    if (success_count > 0)
        return true;
    return false;
}

void Test::Finalize()
{
    JetBase::Finalize();
    cerr << "Accepted: " << fAccepted << endl;
}


bool Test::JetLoop()
{
    unsigned bcount = 0;
    int b1 = 0,b2 = 0, thirties = 0;
    vector<unsigned> flavours;
    for (unsigned i = 0; i < fSortedJets.size(); ++i) {
        if ( i == fJetsPerEvent ) break;
        fJetParts = fSortedJets[i].constituents();

        /* Check the jet flavour if not a generic case */
        if (fDefinition == 1)
            PhysicsFlavor(i);
        else if (fDefinition == 2)
            HadronicFlavor(i);
        else if (fDefinition == 3)
            AlgorithmicFlavor(i);
        else if (fDefinition == 4)
            PhysClusterFlavor(i);
        else if (fDefinition == 5)
            AlgoClusterFlavor(i);

        flavours.push_back(fFlavour);
        if (fFlavour == 5) {
            ++bcount;
            if (bcount == 1)
                b1 = i;
            else if (bcount == 2)
                b2 = i;
        }
        
        if (fSortedJets[i].pt() > 30) {
            ++thirties;
            continue;
        }

        fJetEvent->AddJet(fSortedJets[i].px(),
                          fSortedJets[i].py(),
                          fSortedJets[i].pz(),
                          fSortedJets[i].e(),
                          fJetVars,
                          fWeight,
                          fFlavour);
    }
    if (bcount != 2)
        return false;
    
    PseudoJet t1, t2, t3, t4, t5, t6;
    t1 = fMET + fAuxInputs[fLeptonId];
    if (t1.m() < 60 || t1.m() > 110)
        return false;
    
    vector<PseudoJet> working;
    vector< pair<unsigned,unsigned> > working_idx;
    for (auto i = 0u; i < thirties-1; ++i) {
        if (i == b1 || i == b2) continue;
        if (flavours[i]==0) continue;
        for (auto j = i+1; j < thirties; ++j) {
            if (j == b1 || j == b2) continue;
            if (flavours[j]==0) continue;
            t2 = fSortedJets[i] + fSortedJets[j];
            if (t2.m() > 60 && t2.m() < 110) {
                working.push_back(t2);
                working_idx.push_back( std::make_pair(i,j) );
            }
        }
    }
    if (working.size() == 0)
        return false;
    
    t3 = t1 + fSortedJets[b1];
    t4 = t1 + fSortedJets[b2];
    unsigned success = 0;
    unsigned best;
    for ( unsigned i = 0; i < working.size(); ++i ) {
        t5 = working[i] + fSortedJets[b1];
        t6 = working[i] + fSortedJets[b2];
        if (mass_study(t3.m(),t4.m(),t5.m(),t6.m(),fNoisy)) {
            ++success;
            t2 = working[i];
            best = i;
        }
    }
    if (success == 0)
        return false;
    if (success > 1)
        cerr << "    HOX" << endl;
    
    ++fAccepted;
    
    if (fNoisy) {
        cout << "Lepton W:" << t1.m() << " Jet W:" << t2.m() << endl;
        cout << "Light flavours: " << flavours[working_idx[best].first] << " " << flavours[working_idx[best].second] << endl << endl;
    }
    
    return true;
}