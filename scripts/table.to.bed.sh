#!/bin/bash

input=$1
suffix=.table
input=${input%"$suffix"}
echo $input


tail -n +2 $1 | awk 'OFS="\t"{split($5,a,"|");print $1,$2-1,$2,$3,$4,a[1],a[2]}'  > $input.bed
