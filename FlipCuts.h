// FlipCuts.h
// 1 Dec 2012
// INCLUDE GUARD
#ifndef __FLIPCUTS_H_INCLUDED__
#define __FLIPCUTS_H_INCLUDED__

// Based on FlipEfficiency.h in previous versions of this codes

#include "Pythia.h"                         // Include Pythia headers
#include <fastjet/ClusterSequence.hh>       // fastjet clustering
#include <cmath>                            // for error function
#include <sstream>                          // for string stream
#include <time.h>                           // for random number seed
#include <iostream>                         // for i don't know
#include <iomanip>                          // for setting precision?
#include <fstream>                          // for file in/out
#include <algorithm>                        // for sort

using namespace std;

struct signalregion{
    // this is just a data structure to hold the signal region cuts
    // see e.g. table 2 of SUS-12-017-pas
    unsigned int minJets;
    unsigned int minbJets;
    double minMET;
    double minHT;
    bool plusplus;          // allow same sign + charge leptons
    bool minusminus;        // allow same sign - charge leptons
};

    
/******************************************************************************** 
*   Helper functions that calculate intermediate steps, output, etc.            *
********************************************************************************/

void read_count(vector< pair<string, int> >);
void fill_vector(vector< pair<string, int> > &, string, int);
double get_deltaR(fastjet::PseudoJet, fastjet::PseudoJet);

bool lepton_kinematic_cut(pair<int, fastjet::PseudoJet>);
bool jet_kinematic_cut(pair<int, fastjet::PseudoJet>);
bool lepton_selection_cut(pair<int, fastjet::PseudoJet>);
bool lepton_ID_eff(pair<int, fastjet::PseudoJet>);


// This is the old function. We'll overload the definition and then depreciate
bool lepton_iso_eff(    pair<int, fastjet::PseudoJet>, 
                        vector< pair<int, fastjet::PseudoJet> >);

// This is the revised function. 
bool lepton_iso_eff(    unsigned int,
                        vector< pair<int, fastjet::PseudoJet> >,
                        vector< pair<int, fastjet::PseudoJet> >);
// arguments: lepton array index, lepton array, parton array
// for including leptons into the cone, but not the cone lepton itself
                        
                        
bool b_selection_efficiency(pair<int, fastjet::PseudoJet>);
bool lepton_trig_efficiency(vector< pair<int, fastjet::PseudoJet> >);
bool METefficiency(double, double);
bool HTefficiency(double, double);

bool isLepton(int);

void fill_signalregions(vector<signalregion>&);


vector<pair<int,fastjet::PseudoJet> > apply_cut(
    bool(*)(pair<int,fastjet::PseudoJet>),
    vector<pair<int,fastjet::PseudoJet> >
    );

vector<pair<int,fastjet::PseudoJet> > apply_iso(
    vector<pair<int,fastjet::PseudoJet> >,
    vector<pair<int,fastjet::PseudoJet> >);

bool pTordered(pair<int,fastjet::PseudoJet>, 
    pair<int,fastjet::PseudoJet>);


// END INCLUDE GUARD
#endif __FLIPCUTS_H_INCLUDED__

