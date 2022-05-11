#!/bin/bash

# list of files to be merged is input

samtools merge -b $1 $1.merged.bam
