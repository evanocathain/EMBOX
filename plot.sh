#!/bin/bash

rm plot.gp
echo "set term x11" > plot.gp
np=$1
npminusone=`eval 'echo "$np - 1" | bc'`
colour=`head -1 charges | awk '{if ($1==-1) print "-1"; if ($1==1) print "3"}'`
echo "splot 'positions' every ${np}::0 lt $colour notitle" >> plot.gp
for i in `seq 1 $npminusone`
do
    iplusone=`eval 'echo "$i + 1" | bc'`    
    colour=`head -$i charges | tail -1 | awk '{if ($1==-1) print "-1"; if ($1==1) print "3"}'`
    ## Negative in black, positive in blue ##
    echo "replot 'positions' every ${np}::${i} lt $colour notitle" >> plot.gp
done

gnuplot -persist plot.gp



