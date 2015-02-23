// FlipCommandFileFixer.h
// For modifying Pythia Command Files
// INCLUDE GUARD
#ifndef __FLIPCOMMANDFILEFIXER_H_INCLUDED__
#define __FLIPCOMMANDFILEFIXER_H_INCLUDED__

#include <string>              
#include <sstream>              // for string stream
#include <iostream>             // for i don't know
#include <iomanip>              // for setting precision?
#include <fstream>              // for file in/out
//
#include <algorithm>            // These four are all from
#include <functional>           //  http://stackoverflow.com/
#include <cctype>               //  questions/216823/whats-the-
#include <locale>               //  best-way-to-trim-stdstring
using namespace std;            

bool  FixSpectrum(                    // TRUE if spc file changed successfully
        std::string &templatefile,    // full path of the template spectrum
        std::string &filename,        // full path for output spectrum file
        std::string &blockname,       // block where we're making a replacement
        std::string &blockdivider,    // keyword that divides blocks, e.g. 'BLOCK'
        std::string &lineID,          // identifying
        std::string &newline);
//
// Usage: say you want to fix the mass for particle '1000021' in the
//  spectrum 'myspec.spc' with the value '6.0E+02'. You want to base this
//  spectrum on 'template.spc.' Specify that the divider between blocks is 
//  the string 'BLOCK' and that the specific block you want is 'BLOCK MASS'. 
//  Then input the following:
//      templatefile    = 'template.spc'
//      filename        = 'myspec.spc'
//      blockname       = 'BLOCK MASS'
//      blockdivider    = 'BLOCK'
//      lineID          = '1000021'
//      newline         = '6.0E+02'
//
//      FixCommandFile( 'template.spc', 'myspec.spc', 'BLOCK', 'BLOCK MASS', 
//                      '1000021', '6.0E+02')
//
//  The function will start looking for '1000021' after it hits 'BLOCK
//  MASS' and will continue searching until the next 'BLOCK'.
//  Since the mass lives in BLOCK MASS


bool  FixCommand(                     // TRUE if successful
        std::string &templatefile,    // full path of the template command
        std::string &filename,        // full path for output command file
        std::string &linestart,       // line to change starts with this ...
        std::string &newline);        // ... and is replaced with this line
//
// Usage: analogous to FixSpectrum. The idea is that your program outputs a 
//  spectrum file with a specified name. You can change this name in the 
//  program, but you might forget to update the command file which references
//  this filename. FixCommand can be used to go through a file and search for
//  a line that starts with "linestart" (e.g. "SLHA:file = ") and replace the
//  entire line with something else (e.g. "SLHA:file = different.spc").



// END INCLUDE GUARD
#endif __FLIPCOMMANDFILEFIXER_H_INCLUDED__

