#!/bin/bash

input=$1 ## bam file
out_dir=$2 ## fully specified out file

/home/groups/Spellmandata/heskett/replication.rnaseq/liftover/liftOver $1 \
  /home/groups/Spellmandata/heskett/replication.rnaseq/liftover/hg38ToHg19.over.chain.gz \
  $2 liftover.errors.txt
