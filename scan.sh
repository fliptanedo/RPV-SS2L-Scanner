#!/bin/bash
# HOW TO USE
# arguments are as follows:
# $1 stop   mass starting value
# $2 stop   mass increment
# $3 stop   mass number of steps
# ------------------------------
# $4 gluino mass starting value
# $5 gluino mass increment
# $6 gluino mass number of steps
# ------------------------------
# $7 signal region
# 
for i in $(eval echo {0..$3});
    do
    for j in $(eval echo {0..$6});
        do
        # nice ./RPVgPoint $(( $1 + $2*$i )) $(($4 + $5*$j )) $7 
        echo $(( $1 + $2*$i )) $(($4 + $5*$j )) $7
        ./RPVgPoint $(( $1 + $2*$i )) $(($4 + $5*$j )) $7 
        done
    done 

