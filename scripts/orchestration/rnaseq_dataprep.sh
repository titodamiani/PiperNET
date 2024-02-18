#!/bin/bash

#Cluster proteome
python3 03_cluster_proteome.py

#Parse transcriptome.pep and merge annotations
python3 04_rnaseq_parse_proteome.py

#Shorten and format headers
python3 05_rnaseq_format_headers.py
