#!/usr/bin/env zsh

for i in `ls *.ent`
do
    echo $i
    oname=`basename -s .ent $i`
    p_wat_st -f $i -o $oname.dat
done
