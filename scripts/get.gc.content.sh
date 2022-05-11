#!/bin/bash

srun --mem 30000 -t 200 -c 1 bedtools nuc -fi ../../myron_refs/human_g1k_v37.fasta -bed ../bouhassira.repliseq.dec.20/human_g1k_v37_100kb_s100kb.bed | cut -f 1-3,5 > gc.content.100kb.s100kb.bed


