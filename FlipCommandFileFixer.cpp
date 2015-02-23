/******************************************************************************** 
*   FlipCommandFileFixer.cpp by Flip Tanedo (pt267@cornell.edu)                 *
*   Code for RPVg project, Aug 2012                                             *
*   Contains functions used for modifying spc and cmnd files                    *
********************************************************************************/

#include "FlipCommandFileFixer.h"

bool  FixSpectrum(               // TRUE if spc file changed successfully
        std::string &templatefile,  // full path of the template spectrum
        std::string &filename,      // full path for output spectrum file
        std::string &blockname,     // block where we're making a replacement
        std::string &blockdivider,  // keyword that divides blocks, 'BLOCK'
        std::string &lineID,        // line to change starts with (particle ID)
        std::string &newline){      // line to replace previous

    bool success = false;       // did the replacement work?

    // FILE STREAMS
    ifstream instream;
    instream.open(templatefile.c_str());
    ofstream outstream;
    outstream.open(filename.c_str());
    
    // for internal use in the following loops
    string line;
    string trimmed;
    bool foundblock(false); // found the string blockname
    bool nextblock(false);  // after finding blockname, ran into next block
    
    while (!instream.eof()){
        getline(instream,line);
        
        size_t startpos = line.find_first_not_of(" \t");
        // find the position of the first non-trivial string character
        
        if( string::npos != startpos ){         // if string is nontrivial
            trimmed = line.substr( startpos );  // new string with no whitesp
            
            
            
            // Arrived at the end of the correct block
            if( (trimmed.substr(0,blockdivider.length())==blockdivider) && 
                foundblock ){
                nextblock = true;
                // cout << endl << "FOUND END BLOCK" << endl << endl;
            }
            
            // Found the correct block
            if(trimmed.substr(0,blockname.length())==blockname){
                foundblock = true;
                // cout << endl << "FOUND " << blockname << "!" << endl << endl;
            }
            
            
            // Found correct element within the correct block
            if( (trimmed.substr(0,lineID.length())==lineID)  && 
                (foundblock && !nextblock)
                ){
                    line = newline;
                    // cout << endl << "FOUND IT!" << endl << endl;
                    success = true;
                }
        }
        
        outstream << line << '\n';
                
    }
    
    
    // CLEAN UP FILE STREAMS
    instream.close();
    outstream.close();
    
    return success;
}


bool  FixCommand(                     // TRUE if successful
        std::string &templatefile,    // full path of the template command
        std::string &filename,        // full path for output command file
        std::string &linestart,       // line to change starts with this ...
        std::string &newline){        // ... and is replaced with this line

    bool success = false;       // did the replacement work?

    // FILE STREAMS
    ifstream instream;
    instream.open(templatefile.c_str());
    ofstream outstream;
    outstream.open(filename.c_str());
    
    // for internal use in the following loops
    string line;
    string trimmed;
    // bool foundblock(false); // found the string blockname
    // bool nextblock(false);  // after finding blockname, ran into next block
    
    while (!instream.eof()){
        getline(instream,line);
        
        // If you find a line starting with linestart, replace with newline
        if( line.substr(0,linestart.length()) == linestart){
            line = newline;
            success = true;
        }

        // Put the line in the new command file
        outstream << line << '\n';
                
    }
    
    
    // CLEAN UP FILE STREAMS
    instream.close();
    outstream.close();
    
    return success;
            
}

