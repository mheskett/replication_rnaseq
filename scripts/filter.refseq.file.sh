#!/bin/bash

# header
#chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

tail -n +2 | awk 'OFS="\t"{$1,$3,$4,$11,$10,$2}' 
