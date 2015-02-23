/******************************************************************************** 
*   FlipEfficiencySignal.cpp by Flip Tanedo (pt267@cornell.edu)                 *
*   Code for RPVg project, Aug 2012                                             *
*   Modified for BG, b-tagging Oct 2012                                         *
*   Contains routine for calculating signal efficiency                          *
*                                                                               *
*   What it does:   inputs Pythia object (initialized to a given stop-gluino    *
*                   parameter point) and (a) outputs number of events that      *
*                   pass cuts (b) fills a vector with a list of intermediate    *
*                   counts.                                                     *
*                                                                               *
*   How it works                                                                *
*       1.  Input a Pythia object. This may have itself been initialized        *
*           using a signal parameter point or a background LHE file.            *
*       2.  Initialize a list of intermediate counts. These will be output      *
*           to the counts vector for sanity checks.                             *
*       3.  Generate events. (Loop)                                             *
*       4.  For each event, loop through event particles and keep track of      *
*           leptons, (b) jets, and MET. These objects are stored in vectors     *
*           of particle ID and 4-momentum. Note: b-tagging is treated using     *
*           pythia.process instead of pythia.event since we want to apply the   *
*           CMS efficiencies at parton level.                                   *
*       5.  Go through the stored particle lists and apply efficiencies and     *
*           cuts as prescribed by the CMS paper. Apply signal region cuts       *
*           using the data stored when we filled the signal regions. Note       *
*           that the signal region data is hard-coded into FlipCuts.cpp.        *
*       6.  Apply any jet substructure or additional cuts on the events.        *
*           You'll want to copy and paste a new version of this function when   *
*           playing with different kinds of additional cuts. Don't forget to    *
*           declare these new functions in FlipApplyCuts.h                      *
*       7.  Output number of events that passed all this rigmarole.             * 
*                                                                               *
********************************************************************************/

        /************************************************************************
        * IDENTIFICATION & SELECTION EFFICIENCIES                               *
        * ---------------------------------------                               *
        * Note: selection efficiency = ID efficiency * isolation efficency      *
        * CMS gives selection efficiencies, but since we're recasting to higher *
        * energies we want to account for the effect of boosts on the isolation *
        * efficiency. Thus we're going to kludge this by separating selection   *
        * efficiency into ID and isolation efficiencies. Mike has checked that  *
        * the ID efficiency is basically a constant function of lepton energy.  *
        ************************************************************************/


#include "FlipApplyCuts.h"


int recast(
    Pythia8::Pythia& pythia,                // Pythia object
    vector< pair<string, int> > &counts,    // intermediate data (for checking)
    int iSR,                                // Signal Region #
    int nEvent                              // # events in the Pythia object
    ){
    
    Pythia8::Event& event = pythia.event;       
    Pythia8::Event& process = pythia.process;   
    int nAbort = pythia.mode("Main:timesAllowErrors");
    
    vector<signalregion> signal_region;    // as defined in SUS-12-017 Table 1
    fill_signalregions(signal_region);     // fills data from above paper
    
    // DEBUGGING
    // ofstream debug;
    // debug.open("debug.txt"); // append to end of file
    
    
    /****************************************************************************
    *   SET UP COUNTERS FOR SANITY CHECK COUNTS                                 *
    ****************************************************************************/
    
    int nGenerated  = 0; // # generated events 
    int nKinematic  = 0; // # events that pass kinematic cuts on leptons
    int nLepID      = 0; // # events that pass lepton ID efficiencies
    int nLepIso     = 0; // # events that pass lepton isolation efficiencies
    int nbjetSelect = 0; // # events that pass bJet selection efficiencies
    int nDilepton   = 0; // # events that pass dilepton requirement
    int nDilepTrig  = 0; // # events that pass dilep req & trig efficiency
    int nSS2L       = 0; // # events that pass same sign leptons requirement
    int nJets       = 0; // # events with mininum number of jets
    int nbJets      = 0; // # events with minimum number of tagged b jets
    int nMET        = 0; // # events that pass minimum MET requirement
    int nHT         = 0; // # events that pass minimum HT requirement
    int nCharge     = 0; // # events that pass ++ or -- requirement
    int nPassed     = 0; // # events that passed all cuts
    
    
    /****************************************************************************
    *   GENERATE EVENTS & IMPOSE CUTS                                           *
    ****************************************************************************/
    
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // loop over events
        
        if (!pythia.next()) {                   // if no new event
            if (++iAbort < nAbort) continue;    // if not over abort limit
            cout << " Event generation aborted prematurely, owing to error!\n"; 
            break;
        } // End of 'if no new event'        
        
            
        /************************************************************************
        * SET UP EVENT DATA FOR LATER INSPECTION                                *
        ************************************************************************/
        
        vector< pair<int, fastjet::PseudoJet> > leptons; // generated leptons
        vector< pair<int, fastjet::PseudoJet> > partons; // generated partons
        vector< pair<int, fastjet::PseudoJet> > bpartons;   // b quarks (parton)
        vector< pair<int, fastjet::PseudoJet> > hadrons;    // hadrons in event
        fastjet::PseudoJet METvec (0.0, 0.0, 0.0, 0.0);     // cumulative MET
        double MET (0.0);                                   // MET scalar
        double HT (0.0);                                    // HT scalar
            
        
        /************************************************************************
        * LOOP THROUGH EVENT PARTICLES                                          *
        ************************************************************************/
        
        grabEvent(event, leptons, hadrons);
        grabProcess(process, METvec, partons, bpartons);
        nGenerated++;        
        
        
        /************************************************************************
        * IMPOSE KINEMATIC CUTS AND ID EFFICIENCIES                             *
        ************************************************************************/        
                        
        leptons = apply_cut(lepton_kinematic_cut, leptons);
        if (leptons.size() > 1) nKinematic++; else continue;

        partons = apply_cut(jet_kinematic_cut, partons);
        
        leptons = apply_cut(lepton_ID_eff, leptons);
        if (leptons.size() > 1) nLepID++; else continue;
                
        leptons = apply_iso(leptons, hadrons);
        if (leptons.size() > 1) nLepIso++; else continue;
        
        // Order leptons by pT: do this AFTER isolation since we re-order
        sort (leptons.begin(), leptons.end(), pTordered);

        bpartons = apply_cut(b_selection_efficiency, bpartons);
        if (bpartons.size() > 1) nbjetSelect++; else continue;
        
        if (leptons.size() < 2) continue; else nDilepton++;
        
        if (!lepton_trig_efficiency(leptons)) continue; else nDilepTrig++;
        
        // Same-sign dileptons
        // -------------------
        if (leptons[0].first/abs(leptons[0].first) != 
            leptons[1].first/abs(leptons[1].first)) continue;
        else nSS2L++;
        // Note: assuming that you're only looking at two hardest leptons
        
        
        
        // Signal region cuts: from input
        // ------------------------------
        
        if (partons.size() < signal_region[iSR].minJets) continue;        
        else nJets++;
        
        if (bpartons.size() < signal_region[iSR].minbJets) continue;
        else nbJets++;
        
        MET = METvec.pt(); 
        if (!METefficiency(MET,signal_region[iSR].minMET)) continue;
        else nMET++;
        
        for(unsigned int iPar = 0; iPar < partons.size(); iPar++){
            HT += partons[iPar].second.pt();
        } // end for loop over partons
        // 
        if (!HTefficiency(HT,signal_region[iSR].minHT)) continue;
        else nHT++;
        
        bool minmin = (leptons[0].first > 0) && signal_region[iSR].minusminus;
        bool pluplu = (leptons[0].first < 0) && signal_region[iSR].plusplus;
        
        if (!(minmin || pluplu)) continue;
        else nCharge++;
        
        // Made it this far? YOU PASS
        nPassed++;
        
    } // end for loop, going through Events
    
    
    
    // Fill counts
    // -----------
    fill_vector(counts, "Generated events \t", nGenerated);
    fill_vector(counts, ">1 lep. kin. cuts\t", nKinematic);
    fill_vector(counts, ">1 lep. ID. eff.\t", nLepID);
    fill_vector(counts, ">1 lep. Iso. eff.\t", nLepIso);
    fill_vector(counts, ">1 bjets tagged \t", nbjetSelect);
    fill_vector(counts, "at least two leptons \t", nDilepton);
    fill_vector(counts, "triggered two leptons \t", nDilepTrig);
    fill_vector(counts, "same sign dileptons \t", nSS2L);
    
    // The following cuts depend on the signal region, so we have to
    //  "dynamically" generate their labels
    
    stringstream nJetComment;
    nJetComment << "at least " << signal_region[iSR].minJets << " jets \t";
    fill_vector(counts, nJetComment.str(), nJets);
    
    stringstream nbJetComment;
    nbJetComment << "at least " << signal_region[iSR].minbJets << " b jets \t";
    fill_vector(counts, nbJetComment.str(), nbJets);
    
    stringstream nMETComment;
    nMETComment << "at least " << signal_region[iSR].minMET << " GeV MET \t";
    fill_vector(counts, nMETComment.str(), nMET);
    
    stringstream HTComment;
    HTComment << "at least " << signal_region[iSR].minHT << " GeV HT \t";
    fill_vector(counts, HTComment.str(), nHT);
    
    stringstream nChargeComment;
    if ( signal_region[iSR].minusminus && !signal_region[iSR].plusplus)
        nChargeComment << "only -- leptons \t";
    else if ( !signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "only ++ leptons \t";
    else if ( signal_region[iSR].minusminus && signal_region[iSR].plusplus)
        nChargeComment << "either ++ or -- leptons";
    else nChargeComment << "You fucked up, neither ++ or -- leptons ";
    
    fill_vector(counts, nChargeComment.str(), nCharge);
    
    
    // DEBUGGING
    // debug.close();
        
    return nPassed; 
    
    
} // end void signal_efficiency(...)




// HELPER FUNCTIONS


void grabEvent(Pythia8::Event& event,           // Pythia.event
    vector< pair<int,fastjet::PseudoJet> >& leptons,    // leptons
    vector< pair<int,fastjet::PseudoJet> >& hadrons     // hadrons
    ){
        
    for (int iPart = 0; iPart < event.size(); iPart++){
            
        // Skip things that we can't see
        // -----------------------------
        if (!event[iPart].isFinal()) continue;
        if (!event[iPart].isVisible()) continue;   
        if (abs(event[iPart].eta()) >= 5.0) continue;
                    
        // Save particle 4-momentum and negative contribution to MET
        // ---------------------------------------------------------
        fastjet::PseudoJet momentum(event[iPart].px(),
                                    event[iPart].py(),
                                    event[iPart].pz(),
                                    event[iPart].e());
        
        
        // Keep track of leptons
        // ---------------------
        if ( isLepton( event[iPart].id() ) ) {
        
            pair<int,fastjet::PseudoJet> cur_lept;
            cur_lept.first = event[iPart].id();     // PDG code
            cur_lept.second = momentum;             // 4-vector
            leptons.push_back(cur_lept);         // put in list
            continue;
            
        } // End "if this is an identfiable lepton"  
        
        // Leftover objects are non-leptonic, i.e. hadrons
        pair<int,fastjet::PseudoJet> cur_hadron;
        cur_hadron.first = event[iPart].id();     // PDG code
        cur_hadron.second = momentum;             // 4-vector
        hadrons.push_back(cur_hadron);
        
        } // End loop through event particles
        
    }
    

void grabProcess(Pythia8::Event& process,   // Pythia.process
    fastjet::PseudoJet& METvec,             // METvec
    vector< pair<int,fastjet::PseudoJet> >& partons,   // partons
    vector< pair<int,fastjet::PseudoJet> >& bpartons   // bpartons
    ){
        
    for (int iPart = 0; iPart < process.size(); iPart++){    
        if (!process[iPart].isFinal()) continue;
        if (!process[iPart].isVisible()) continue;   
        if (abs(process[iPart].eta()) >= 5.0) continue;
        
        fastjet::PseudoJet momentum(process[iPart].px(),
                                    process[iPart].py(),
                                    process[iPart].pz(),
                                    process[iPart].e());
                                    
        METvec -= momentum; // generator level MET for recast
        
        if ( isLepton( process[iPart].id() ) ) continue; // ignore e, mu
                                    
        // Anything left is a parton
        // -------------------------    
        pair<int,fastjet::PseudoJet> cur_parton;
        cur_parton.first = process[iPart].id();
        cur_parton.second = momentum;
        partons.push_back(cur_parton); 
    
        // b quarks
        // --------
        if (!(abs(process[iPart].id())==5)) continue;     // only bjets
        pair<int,fastjet::PseudoJet> cur_bparton;
        cur_bparton.first = process[iPart].id();
        cur_bparton.second = momentum;
        bpartons.push_back(cur_bparton); 
        
    } // End loop through process particles
        
} // end grabProcess


