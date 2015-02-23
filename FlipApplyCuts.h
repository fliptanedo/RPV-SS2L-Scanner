// FlipApplyCuts.h
// 1 Dec 2012
// INCLUDE GUARD
#ifndef __FLIPAPPLYCUTS_H_INCLUDED__
#define __FLIPAPPLYCUTS_H_INCLUDED__

// Based on FlipEfficiency.h in previous versions of this codes

#include "Pythia.h"                         // Include Pythia headers
#include <fastjet/ClusterSequence.hh>       // fastjet clustering
#include <cmath>                            // for error function
#include <sstream>                          // for string stream
#include <time.h>                           // for random number seed
#include <iostream>                         // for I don't know
#include <iomanip>                          // for setting precision?
#include <fstream>                          // for file in/out
#include "FlipCuts.h"                       // for cut/efficiency tools
using namespace std;

int recast(Pythia8::Pythia&, vector< pair<string, int> >&, int, int);
    // This is our main workhorse, it's defined in FlipApplyCuts.cpp
    // Inputs: pythia object, count vector, signal region index, # event
    // Output: number of events that pass the cuts

// Eventually we'll want to have different kinds of functions
// E.g. for doing substructure, etc.


// HELPER FUNCTIONS

void grabEvent(Pythia8::Event&,                 // Pythia.event
    vector< pair<int,fastjet::PseudoJet> >&,    // leptons
    vector< pair<int,fastjet::PseudoJet> >&     // hadrons
    );

void grabProcess(Pythia8::Event&,   // Pythia.process
    fastjet::PseudoJet&,            // METvec
    vector< pair<int,fastjet::PseudoJet> >&,   // partons
    vector< pair<int,fastjet::PseudoJet> >&    // bpartons
    );



// END INCLUDE GUARD
#endif __FLIPAPPLYCUTS_H_INCLUDED__

