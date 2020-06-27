#!/bin/bash

Rscript R_scripts/rocAnalysis_for_Ritu.R -i Data/Fasta/mod_scerevisiae.fasta -o Scer/Example_run_for_Ritu/ -d 0 -s 1000 -a 50 -t 2 -n 1 --est_csp --est_phi --est_hyp &> example_run_for_Ritu.Rout