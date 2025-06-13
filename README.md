<p align="center">
<img src=https://github.com/user-attachments/assets/fb39c9dc-4803-4684-b65a-935b9202c5a9 width=300 height=140 align=center>
</p>

# Alag
Alag is a comparative and population genomics tool that performs in-depth analyses of the action of constraint and adaptation on the **protein-coding** elements within the genomes

## Prerequisites

Alag is exclusively written in R, with occasional internal calls to external tools (vcftools, blastn and MUSCLE) that are to be preinstalled. The R libraries required for running *Alag* are listed in **key.R**.

## Structure

The tool's centre point is **key.R**. This script is the import point for libraries and the species-specific input files that will be used later for the analyses. **main.R** contains the call points for scripts (in the __includes__ folder) that process the input files, which are passed on sequentially to the consecutive functions.

## How to run

The input points for the species-specific files are documented in **key.R**. These can be replaced with input files for other species. Example input files are listed in the folder __input_files__ . The resulting output files are listed in __example_outputs__ folder. Scripts used for performing analyses are listed in __analysis__ folder

