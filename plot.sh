#!/bin/bash

rm plot.gp
echo "set term x11" > plot.gp
np=$1
npminusone=`eval 'echo "$np - 1" | bc'`
echo "splot 'positions' every ${np}::0 notitle" >> plot.gp
for i in `seq 1 $npminusone`
do
    echo $i
    echo "replot 'positions' every ${np}::${i} notitle" >> plot.gp
done

gnuplot -persist plot.gp

#splot './positions' every np::0 lt 0, './positions' every np::1 lt 0, './positions' every np::2 lt 0, './positions' every np::3 lt 0, './positions' every np::4 lt 0, './positions' every np::5 lt 0, './positions' every np::6 lt 0, './positions' every np::7 lt 0, './positions' every np::8 lt 0, './positions' every np::9 lt 0

