SS2L search for gluino decays into light RPV stops
RECAST code by Flip Tanedo (validating code by Mike Saelim)
with Mike Saelim, Maxim Perelstein, and Josh Berger

13 Jan 2013
pt267@cornell.edu


WHAT DOES IT DO: This code calculates the signal efficiency for SS2L production
at the LHC in a simplified model with stops and gluinos and an RPV coupling
W ~ UDD. The output is stored in an output file (output.dat) and the code
is written to make it easy to scan over parameters using a bash script.

This code has been tested to work with Pythia 8.165.

FILES:
------
Makefile:   Makefile
Drivers:    RPVgPoint.cc
Templates:  TEMPLATE.spc
            TEMPLATE.cmnd
Auxiliary:  FlipCommandFileFixer.cpp/h
            FlipCuts.cpp/h
            FlipApplyCuts.cpp/h
Output:     output.dat
Temporary:  TEMP.spc
            CommandRun.cmnd
            spcRun.spc

(temporary files are created and deleted each run)


USAGE:
------

1. Go to the Makefile and modify the first two non-commented lines
    These correspond to the path to your PYTHIA and FASTJET folders respectively
    Everything else in this file should work out-of-the-box.
    You might want to change the compiler/flags. Whatever floats your boat, dude.
   
    
2. In the terminal, compile this program using the command 'make'
    If there are errors, then it's either your fault or my fault.
   
    
3. If successful, make will output instructions for how to use the compiled code.
    There are 'template' command and spectrum files which the program runs. when
    no inputs are specified. Please check RPVgPoint.cc for defaults. 
    
    
    
    (Default template files CAPITALIZED filenames)
    
    If you want to modify something in the Pythia run, go ahead and modify the
    command file template (or even better, create a new one). For example, the
    current code assumes that we're forcing leptons to decay leptonically. 
    (The code compensates for this with an overall prefactor to the efficiency
    that is sent to the output file. Note that the efficiency that is listed
    upon running does not yet have this prefactor.) You might want to change this
    in the template command file.
    
4. Anyway, go ahead and run the program following the sample call from the
    instructions listed in the makefile. For example, you can run with default
    values without any options:
    
        ./RPVgPoint
    
    You can also put in a bunch of options:

        ./RPVgPoint 300 800 8 CmndShort.cmnd output.dat template.spc
    
    These correspond to (in order):
        stop mass
        gluino mass
        signal region (for comparison with CMS SUS-12-017)
        command file template (CmndShort.cmnd only runs 100 events)
        output file
        spectrum file template

    In other words:
	./RPVgPoint [mstop] [mglu] [SigReg] [cmnd] [output] [spc]
    
    There are default values for each of these, but if you specify any options,
    you have to make sure all possible options before are also specified. In 
    other words, if you only wanted to change the signal region, you would still
    have to specify the stop and gluino masses, since, for example,
    
        ./RPVgPoint 8
        
    would be interpreted as setting the stop mass to 8.
    
5. Scanning with a batch script: this was the raison d'etre for this code. 
    This is straightforward since you can just scan over the program options.
    The script scan.sh automates the scan once RPVgPoint is compiled. Make sure
    that scan.sh has the appropriate permissions to be executed (use the chmod
    command). The script takes 7 arguments:
	$1 stop   mass starting value
	$2 stop   mass increment
	$3 stop   mass number of steps
	$4 gluino mass starting value
	$5 gluino mass increment
	$6 gluino mass number of steps
	$7 signal region
    So, for example, one can run:
	./scan.sh 200 10 3 1200 10 3 8
    You have to modify scan.sh directly if you want to change the other options,
    e.g. if you want to use different template cmnd or spc files.
    
    
    
MORE DETAILS ON HOW THE CODE WORKS

The driver program (e.g. PartonLevel.cc) creates a Pythia object and collects the
relevant run data. It then passes the object and run data into a function ApplyCuts
which calls several helper functions to apply kinematic cuts and selection 
efficiencies to the Pythia events. The output of ApplyCuts is a number of events
that passes the cuts. 

There is also a vector "counts" which stores the number of events after each cut.
This is mainly for debugging to check if certain cuts are misbehaving. This vector
is filled in the section "Fill counts" of FlipApplyCuts.cpp.

Note that the RPVgPoint.cc code has a default command file that asks Pythia to do
the full hadronic event, but most of the code uses only the parton-level data found
in pythia.process (this is an "event" type object). The reason for this is that 
CMS gives lepton selection efficiencies, which is equal to the product of the ID
efficiency and the isolation efficiency. Unfortunately, the isolation efficiency is
a function of the boost of the decaying heavy object, and so we should calculate this
ourselves. (Mike has checked that the ID efficiency function is basically constant.)
In calculating the lepton isolation, it is necessary to consider the hadronic
activity that might fall into a lepton's isolation cone and so it makes sense to use
the full hadronic event rather than the parton-level event with 'pencil jets'.



REMARKS AND THINGS TO DO

1. Note that the "same sign dilepton" check only looks at the two hardest leptons.
    In principle one might want to code up a more comprehensive function to check this.
    
    
    
Good scanning,
Flip, Sept 2012 (rev Jan 2013)
    