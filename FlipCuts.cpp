/******************************************************************************** 
*   FlipEfficiency.cpp by Flip Tanedo (pt267@cornell.edu)                       *
*   Code for RPVg project, Aug 2012                                             *
*   Contains functions used when calculating signal efficiency                  *
********************************************************************************/

#include "FlipCuts.h"



void read_count(vector< pair<string, int> > count){
    // outputs the contents of count to screen
    
    cout << endl;
    for(unsigned int i=0; i < count.size(); i++)
        cout << count[i].first << ":\t" << count[i].second << endl;
} // end void read_count(...)



void fill_vector(vector< pair<string, int> > &count, string line, int num){
    // adds an element to a vector of particles
    
    pair<string, int> new_item(line, num);
    count.push_back(new_item);
} // end void fill_vector(...)



double get_deltaR(fastjet::PseudoJet vec1, fastjet::PseudoJet vec2){
    // outputs the Delta_R between two four-momenta (pseudoJets)

    double phi1 = vec1.phi();
    double eta1 = vec1.eta();
    double phi2 = vec2.phi();
    double eta2 = vec2.eta();

    double Rsq = pow(phi2-phi1,2) + pow(eta2-eta1,2);
    return sqrt(Rsq);
} // end get_deltaR



bool lepton_kinematic_cut(pair<int, fastjet::PseudoJet> lepton){
    // returns true if a lepton passes the kinematic cuts
    
    // LEPTON KINEMATIC CUT PARAMETERS
    double electron_pT  = 20.0;
    double muon_pT      = 20.0;
    double lepton_eta   = 2.4;
    double eta_bar   = 1.442;
    double eta_end   = 1.566;
    
    bool passes = false;
    int id = lepton.first;
    fastjet::PseudoJet momentum = lepton.second;
    
    bool pass_pT =  ((abs(id) == 11) && (momentum.pt() >= electron_pT)) ||
                    ((abs(id) == 13) && (momentum.pt() >= muon_pT));
                    
    bool pass_eta = ((abs(id) == 13) && (abs(momentum.eta()) < lepton_eta)) ||
                    ((abs(id) == 11) && (abs(momentum.eta()) < eta_bar)) ||
                    ((abs(id) == 11) && (abs(momentum.eta()) > eta_end)
                                     && (abs(momentum.eta()) < lepton_eta));
                                     
    if (pass_pT && pass_eta) passes = true;
    
    return passes;
        
} // end lepton_kinematic_cut



bool jet_kinematic_cut(pair<int, fastjet::PseudoJet> jet){
    // returns true if a jet passes the kinematic cuts
    // note that in SUS-12-017 the jet and bjet kin cuts are the same
    //  so I haven't written a separate bjet_kinematic_cut function
    
    // JET KINEMATIC CUT PARAMETERS
    double jet_pT   = 40.0;
    double jet_eta  = 2.4;
    
    bool passes = false;
    // int id = jet.first; 
    fastjet::PseudoJet momentum = jet.second;
    
    bool pass_pT    = (momentum.pt() >= jet_pT);                    
    bool pass_eta   = (abs(momentum.eta()) < jet_eta);
                                     
    if (pass_pT && pass_eta) passes = true;    
    
    return passes;
        
} // end jet_kinematic_cut



bool lepton_selection_cut(pair<int, fastjet::PseudoJet> lepton){
    // Selection efficiency for leptons
    // note: no longer used in favor of separate ID and iso efficiencies
    
    bool passes = false;
//    int random = rand() % 1001; // random number from 0 to 1000
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    
    // LEPTON EFFICIENCY PARAMETERS
    // Parameterization in eq. 1 of SUS-12-017-pas
    double einf_e   = 0.58;
    double einf_mu  = 0.66;
    double e20_e    = 0.22;
    double e20_mu   = 0.47;
    double sig_e    = 12;
    double sig_mu   = 26;
    double pt = lepton.second.pt();
    
    double einf = 0;
    double e20 = 0;
    double sigma = 1;
    
    if (abs(lepton.first) == 11){
        einf = einf_e;
        e20 = e20_e;
        sigma = sig_e;
    }
    
    if (abs(lepton.first) == 13){
        einf = einf_mu;
        e20 = e20_mu;
        sigma = sig_mu;
    }
    
    double first_term = einf * erf((pt - 20)/sigma);
    double second_term = e20 * erfc((pt - 20)/sigma);
    double efficiency = first_term + second_term;
    
    
    // if (random < 1000*efficiency ) passes = true;
    if (random < efficiency ) passes = true;

    return passes;
    
} // end lepton_selection_cut



bool lepton_ID_eff(pair<int, fastjet::PseudoJet> lepton){
    // Lepton ID efficiency
    
    bool passes = false;
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    
    // LEPTON EFFICIENCY PARAMETERS

    double IDefficiency = 0.0;     
    
    if (abs(lepton.first) == 11) IDefficiency = 0.76;   // electron    
    if (abs(lepton.first) == 13) IDefficiency = 0.86;   // muon
    
    if (random < IDefficiency) passes = true;

    return passes;
    
} // end lepton_ID_eff



bool lepton_iso_eff(    pair<int, fastjet::PseudoJet> lepton, 
                        vector< pair<int, fastjet::PseudoJet> > partons){
    // Lepton isolation efficiency
    
    bool passes = false;
    double cone_pT = 0;
    double lepton_dR  = 0.3; // lepton delta R
    double Iiso       = 0.15;
    
    for (unsigned int iJet = 0; iJet < partons.size(); iJet++) {
        if (get_deltaR(lepton.second, partons[iJet].second) < lepton_dR)
            cone_pT += partons[iJet].second.pt();
    } // end loop over parton
    
    if (cone_pT < Iiso*lepton.second.pt() ) passes = true;
    return passes;

    
                            
} // end lepton_iso_eff



bool lepton_iso_eff(    unsigned int seedLepton,
                        vector< pair<int, fastjet::PseudoJet> > leptons, 
                        vector< pair<int, fastjet::PseudoJet> > partons){
    // Lepton isolation efficiency
    
    bool passes = false;
    double cone_pT = 0;
    double lepton_dR  = 0.3; // lepton delta R
    double Iiso       = 0.15;
    
    pair<int, fastjet::PseudoJet> lepton = leptons[seedLepton];
    
    // Fill cone_pT with cone partons
    for (unsigned int iJet = 0; iJet < partons.size(); iJet++) {
        if (get_deltaR(lepton.second, partons[iJet].second) < lepton_dR)
            cone_pT += partons[iJet].second.pt();
    } // end loop over partons
    
    // Fill cone_pT with cone leptons, don't count seed lepton
    for (unsigned int iLep = 0; iLep < leptons.size(); iLep++) {
        if (get_deltaR(lepton.second, leptons[iLep].second) < lepton_dR){
            if(iLep != seedLepton)
                cone_pT += leptons[iLep].second.pt();
        } // end if lepton is in the cone
    } // end loop over leptons
    
    if (cone_pT < Iiso*lepton.second.pt() ) passes = true;
    return passes;
                            
} // end lepton_iso_eff (3 arguments)



bool b_selection_efficiency(pair<int, fastjet::PseudoJet> bjet){
    // based on efficiencies, randomly determines if
    // a generated bjet is successfully tagged
    
    bool passes = false;
    // int random = rand() % 1001; // random number from 0 to 1000
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    
    
    double pt = bjet.second.pt();
    double efficiency = .65;
    
    // parameterization form SUSY-12-917-pas
    if( (pt > 90) && (pt < 170)) efficiency = 0.65;
    else if (pt <= 90) efficiency = .65 - (90 - pt) * 0.0038;
    else if (pt >= 170 ) efficiency = .65 - (pt - 170) * 0.0007;
    
    if (pt < 40) efficiency = 0; // cut on bjet
    
    // if (random < efficiency*1000) passes = true;
    if (random < efficiency) passes = true;
    
    return passes;
} // end tag_b



bool lepton_trig_efficiency(vector< pair<int, fastjet::PseudoJet> > leptons){
    // Gives probability that a dilepton pair is triggered upon
    // Should also require one lepton with pT > 17, other with pT > 8
    //  but this is already automatically satisfied by lepton kinematic cuts
    // Make sure you sort leptons by decreasing pT so you're testing the
    //  two hardest leptons. See FlipApplyCuts.cpp.
    
    bool passes = false;
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    double eff_ee = 0.95;
    double eff_emu = 0.92;
    double eff_mumu = 0.88;
    
    // if (random < eff_emu) return true;
    // else return false;
    
    if (leptons.size()<2) return passes;   // need at least 2 leptons
    
    if ((abs(leptons[0].first) == 11) && (abs(leptons[1].first) == 11))
        if (random < eff_ee) passes = true;
    if ((abs(leptons[0].first) == 11) && (abs(leptons[1].first) == 13))
        if (random < eff_emu) passes = true;
    if ((abs(leptons[0].first) == 13) && (abs(leptons[1].first) == 11))
        if (random < eff_emu) passes = true;
    if ((abs(leptons[0].first) == 13) && (abs(leptons[1].first) == 13))
        if (random < eff_mumu) passes = true;
    //     
    return passes;
            
    
    // // Minimum trigger pT cuts
    // // Note that these are automatically satisfied by the kinematic cuts
    // double big_pT = 17.0;  // GeV, at least one lepton minimum pT
    // double small_pT = 8.0; // GeV, other lepton minimum pT
    // 
    // if (leptons.size()!=2) return passes;   // exactly two leptons, check
    // 
    // bool bigCut = false;                    // one lepton with pT > upperpT
    // bool smallCut = false;                  // one lepton with pT > lowerpT
    // 
    // // Check that both lepton pTs pass the small cut
    // if (    (leptons[0].second.pt() > small_pT) &&
    //         (leptons[1].second.pt() > small_pT) ) smallCut = true;
    // 
    // // Check that one of the lepton pTs passes the big cut
    // for(unsigned int iLep=0; i< leptons.size(); i++){
    //     if (leptons[iLep].second.pt() > big_pT) bigCut = true;
    // }

    
    
} // end lepton_trig_efficiency



bool METefficiency(double MET, double minMET){
    // Converts between parton-level MET and hadronic MET
    // by including effect of 'turn on curves'
    // from 1205.3933
    
    bool passes = false;
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    
    double x = MET;
    double x12 = 0;
    double sig = 0;
    
    if (minMET >= 120){
        x12 = 123;
        sig = 37;
    }
    else if (minMET >= 50){
        x12 = 43;
        sig = 39;
    }
    else if (minMET >= 30){
        x12 = 13;
        sig = 44;
    }
    else if (minMET == 0); // do nothing, see below
    else cout << endl << "ERROR: METefficiency" << endl;
    
    double efficiency = 0.5*(erf((x-x12)/sig) + 1);
    if (random < efficiency) passes = true;
    
    if (minMET == 0) passes = true;
    
    return passes;
    
} // end METefficiency



bool HTefficiency(double HT, double minHT){
    // Converts between parton-level HT and hadronic HT
    // by including effect of 'turn on curves'
    // from 1205.3933
    
    bool passes = false;
    double random = (double)rand()/(double)RAND_MAX; // random from 0 to 1
    
    double x = HT;
    double x12 = 0;
    double sig = 0;
    
    if (minHT >= 320){
        x12 = 188;
        sig = 88;
    }
    else if (minHT >= 200){
        x12 = 308;
        sig = 102;
    }
    else if (minHT == 80); // do nothing, see below
    else if (minHT == 0);  // equivalent to above cut
    else cout << endl << "ERROR: HTefficiency" << endl;
    
    double efficiency = 0.5*(erf((x-x12)/sig) + 1);
    if (random < efficiency) passes = true;
    
    if ( (minHT == 80) || (minHT == 0) ) passes = true;
    // minimum pT cuts on jet selection is 40 GeV
    // so a min HT of 80 trivially passes cuts
    
    return passes;
} // end HTefficiency


bool isLepton(int pid){
    bool result = false;
    if ( abs(pid) == 11)
        result = true;
    if ( abs(pid) == 13)
        result = true;
        
    return result;
} // end isLepton



vector<pair<int,fastjet::PseudoJet> > apply_cut(
    bool(*pass)(pair<int,fastjet::PseudoJet>),
    vector<pair<int,fastjet::PseudoJet> > particles){
    // 
    
    vector<pair<int,fastjet::PseudoJet> > output;
    for(unsigned int iPar = 0; iPar < particles.size(); iPar++){
        if (pass(particles[iPar])) output.push_back(particles[iPar]);
    } // end for loop over leptons
    
    return output;
}



vector<pair<int,fastjet::PseudoJet> > apply_iso(
    vector<pair<int,fastjet::PseudoJet> > leptons,
    vector<pair<int,fastjet::PseudoJet> > hadrons){
    //
    
    vector< pair<int, fastjet::PseudoJet> > templeptons;
    
    for(unsigned int iLep = 0; iLep < leptons.size(); iLep++){
        if (lepton_iso_eff(iLep, leptons, hadrons)) 
            templeptons.push_back(leptons[iLep]);
    } // end for loop over leptons
    
    return templeptons;   
}



bool pTordered(pair<int,fastjet::PseudoJet> part1, 
    pair<int,fastjet::PseudoJet> part2){
    //
    
    return ( part1.second.pt() >  part2.second.pt());
}
    



void fill_signalregions(vector<signalregion>& signal_region){
    // Fills signal regions with data from PAS SUS-12-029, table 2
    // Input: empty "signalregion" vector
    // Note: see Mike's code for a more efficient way of doing this
    
    signal_region.push_back(signalregion());
    signal_region[0].minJets    = 2;
    signal_region[0].minbJets   = 2;
    signal_region[0].minMET     = 0.0;
    signal_region[0].minHT      = 80.0;
    signal_region[0].plusplus   = true;
    signal_region[0].minusminus = true;
    
    
    signal_region.push_back(signalregion());
    signal_region[1].minJets    = 2;
    signal_region[1].minbJets   = 2;
    signal_region[1].minMET     = 30.0;
    signal_region[1].minHT      = 80.0;
    signal_region[1].plusplus   = true;
    signal_region[1].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[2].minJets    = 2;
    signal_region[2].minbJets   = 2;
    signal_region[2].minMET     = 30.0;
    signal_region[2].minHT      = 80.0;
    signal_region[2].plusplus   = true;
    signal_region[2].minusminus = false;
    // SR2: same as SR1, but only ++ lepton pairs, no --
    // will implement this elsewhere
    
    signal_region.push_back(signalregion());
    signal_region[3].minJets    = 4;
    signal_region[3].minbJets   = 2;
    signal_region[3].minMET     = 120.0;
    signal_region[3].minHT      = 200.0;
    signal_region[3].plusplus   = true;
    signal_region[3].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[4].minJets    = 4;
    signal_region[4].minbJets   = 2;
    signal_region[4].minMET     = 50.0;
    signal_region[4].minHT      = 200.0;
    signal_region[4].plusplus   = true;
    signal_region[4].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[5].minJets    = 4;
    signal_region[5].minbJets   = 2;
    signal_region[5].minMET     = 50.0;
    signal_region[5].minHT      = 320.0;
    signal_region[5].plusplus   = true;
    signal_region[5].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[6].minJets    = 4;
    signal_region[6].minbJets   = 2;
    signal_region[6].minMET     = 120.0;
    signal_region[6].minHT      = 320.0;
    signal_region[6].plusplus   = true;
    signal_region[6].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[7].minJets    = 3;
    signal_region[7].minbJets   = 3;
    signal_region[7].minMET     = 50.0;
    signal_region[7].minHT      = 200.0;
    signal_region[7].plusplus   = true;
    signal_region[7].minusminus = true;
    
    signal_region.push_back(signalregion());
    signal_region[8].minJets    = 4;
    signal_region[8].minbJets   = 2;
    signal_region[8].minMET     = 0.0;
    signal_region[8].minHT      = 320.0;
    signal_region[8].plusplus   = true;
    signal_region[8].minusminus = true;
}

