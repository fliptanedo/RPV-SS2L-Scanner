/******************************************************************************** 
*   PartonLevel.cc by Flip Tanedo (pt267@cornell.edu)                           *
*   Scan of RPV stop/gluino parameter space to determine the SS2L reach         *
*   01 Dec 2012                                                                 *
*                                                                               *
*   This file contains the 'meat' for doing signal and background runs, however *
*   the default configuration is for a signal run. I suggest duplicating and    *
*   editing if you want to run a batch background run. For example, you'll want *
*   to edit the default argc input values to something more appropriate for a   *
*   Standard Model background run.                                              *
*                                                                               *
********************************************************************************/



#include "FlipCommandFileFixer.h"   // to update command file
#include "FlipCuts.h"               // all of my functions
#include "FlipApplyCuts.h"          // all of my functions
#include "Pythia.h"                 // Include Pythia headers
#include <vector>                   // for vectors
#include <sstream>                  // for string stream
#include <fstream>                  // for file in/out
// 
using namespace std;

int main(int argc, char *argv[]) { 

    /****************************************************************************
    * DEFINE AND INITIALIZE PARAMETERS                                          *
    ****************************************************************************/

    // INITIALIZE 
    // ----------
    srand((unsigned)time(0));               // Initialize random numbers
    string outfile = "output.dat";          // Output filename
    vector< pair<string, int> > counts;     // counts @ each cut w/ descriptions
    vector<string> tempfiles;               // Intermediate files to be deleted


    // A bunch of definitions for setting the stop and gluon masses
    // -------------------------------------------------------------
    string mgluino      = "800";                // default gluino mass        
    string mstop        = "300";                // default stop mass
    string cmndtemp     = "TEMPLATE.cmnd";      // default cmnd file template
    string cmndrun      = "CommandRun.cmnd";    // output cmnd file for run
    string spctemp      = "TEMPLATE.spc";       // default spc template
    string spcint       = "TEMP.spc";           // intermediate spc file
    string spcRun       = "spcRun.spc";         // output spc file for run
    //
    string cmndbg       = "TEMPLATEBG.cmnd";    // command file for BG run
    string input_lhe    = "BGevents.lhe";       // input LHE file
    //
    string cmndspc      = "SLHA:file = ";       // line to change in cmnd file
    string cmndspcnew   = cmndspc + spcRun;     // ... replace with this
    //
    string blockmass    = "BLOCK MASS";         // BLOCK MASS tag
    string blockdiv     = "BLOCK";              // BLOCK divider tag
    string gluinoID     = "1000021";            // Gluino PDG code
    string stopID       = "1000006";            // Stop PDG code

    // Identify which intermediate files should be deleted
    // ---------------------------------------------------
    
    tempfiles.push_back(cmndrun);
    tempfiles.push_back(spcint);
    tempfiles.push_back(spcRun);
    
    
    
    // Other definitions for the run
    // -----------------------------
    int iSR = 8;        // Signal region #, defined in SUS-12-017
    
    
    // TAKE IN EXTERNAL VALUES
    // -----------------------
    // You'll want to change these if you're doing a batch run of 
    // background events.
    //
    if (argc > 1)  mstop    = argv[1];       // stop mass
    if (argc > 2)  mgluino  = argv[2];       // gluino mass
    if (argc > 3)  iSR      = atoi(argv[3]); // signal region
    if (argc > 4)  cmndtemp = argv[4];       // template command file
    if (argc > 5)  outfile  = argv[5];       // output filename
    if (argc > 6)  spctemp  = argv[6];       // template spectrum file



    // OUTPUT FILE STREAM
    // ------------------

    ofstream outstream;
    outstream.open(outfile.c_str(), ios::app); // append to end of file
    outstream.precision(6); 
    outstream.setf(ios::fixed);
    outstream.setf(ios::showpoint);



    /****************************************************************************
    * UPDATE SPECTRUM ACCORDING TO PARAMETER SPACE POINT                        *
    * --------------------------------------------------                        *
    * This part of the code uses commands from FlipCommandFileFixer to generate *
    * new command and spectrum files with the desired stop and gluino masses,   *
    * as defined above.                                                         *
    *                                                                           *
    ****************************************************************************/

    // Lines for updating the spectrum
    // --------------------------------
    string gluinoNew    = "   1000021   " + mgluino;    // line replacement
    string stopNew      = "   1000006   " + mstop;      // line replacement
    
    // Make sure command file is using the same spc file that we're creating
    // ---------------------------------------------------------------------
    if(!FixCommand(cmndtemp, cmndrun, cmndspc, cmndspcnew)) 
        cout << endl << " ERROR in FixCommand, setting " << cmndspcnew << endl;

    // Using template spectrum, set gluino mass. Save to intermediate spectrum
    // -----------------------------------------------------------------------
    if(!FixSpectrum(spctemp, spcint, blockmass, blockdiv, gluinoID, gluinoNew))
        cout << endl << "ERROR: FixSpectrum, setting mass " << gluinoNew << endl;
    
    // Using intermediate spectrum, set stop mass. Save to final spectrum
    // -------------------------------------------------------------------
    if(!FixSpectrum(spcint, spcRun, blockmass, blockdiv, stopID, stopNew))
        cout << endl << "ERROR: FixSpectrum, setting mass " << stopNew << endl;



    /****************************************************************************
    * DEFINE AND INITIALIZE PYTHIA OBJECT FOR SHOWERING                         *
    * -------------------------------------------------                         *
    * Here we create the Pythia object and initialize according to whether we   *
    * are calculating signal or background.                                     *
    *                                                                           *
    * SIGNAL:   input cmndrun (generated cmnd file which inputs generated       *
    *           spectrum file) and initialize.                                  *
    * BCKGRND:  Do not use SUSY info above. Instead, initialize on an LHE file  *
    *           generated by MadGraph.                                          *
    *                                                                           *
    * Pick either SIGNAL or BACKGROUND, but don't initialize both on the same   *
    * Pythia object! Instead, define two objects:                               *
    *   Pythia8::Pythia pythiasignal;                                           *
    *   Pythia8::Pythia pythiabackground;                                       *
    *                                                                           *
    ****************************************************************************/

    // SIGNAL INITIALIZATION
    // ---------------------
    Pythia8::Pythia pythia;                     // Declare Pythia object
    pythia.readFile(cmndrun);                   // Read in command file

    int nEvent = pythia.mode("Main:numberOfEvents");
    pythia.init();



    /****************************************************************************
    *   THIS PART DOES THE CALCULATION                                          *
    *****************************************************************************/

    outstream << mstop << "\t" << mgluino << "\t" << iSR << "\t" 
        << double(recast(pythia, counts, iSR, nEvent)) * .10608  
        << "\t" << nEvent << endl;
        // 
        // When calculating efficiency, don't forget to include a factor of
        // 0.10608 = 0.3257^2 from W decays forced to go to leptons (for stats)
        // Why? Because we assume you forced the W to decay leptonically in
        // the command file.
        //
        // NOTE: technically, instead of nEvent, one should use the data from
        //  the counts vector since this 'total generated events' number 
        //  accounts for any aborted events. It shouldn't be a big difference
        //  since we abort the entire run if too many events fail.

    // // IF YOU WANT VERBOSE SCREEN OUTPUT:
    // cout << "STOP: " << mstop << endl;
    // cout << "GLUINO: " << mgluino << endl;
    // cout << "Signal Region " << iSR << endl; 
    read_count(counts); // gives intermediate steps
    cout << endl;
    // cout << endl << endl;   



    /****************************************************************************
    *   CLEAN UP: close filestreams, etc.                                       *
    *****************************************************************************/
    
    outstream.close();


    // REMOVE TEMPORARY FILES
    // ----------------------        
    for(unsigned int iStr = 0; iStr < tempfiles.size(); iStr++){
        if( remove(tempfiles[iStr].c_str()) !=0 )
            cout << "ERROR deleting " << tempfiles[iStr] << endl;
        else
            cout << tempfiles[iStr] << " deleted correctly." << endl;
    } // end loop over temporary files
    
    
    
    return 0;
        
}
    
