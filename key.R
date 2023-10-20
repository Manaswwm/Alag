## centre-point for Alag local code
## this package will replicate all the Alag functions on local machines

# ================    INPUT FILES   ==========================================================================
#gff file for ingroup species; note: structure is quite important for that one! see function gff_extraction.R
ingroup_gff_path="/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/gff_aratha/Arabidopsis_thaliana.TAIR10.55.gff3"

#fasta files with transcript sequences
ingroup_transcript_cds_path =  "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/ingroup_aratha/aratha_transcripts_proteincoding.fa"
outgroup_transcript_cds_path= "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/outgroup_alyra/alyra_transcripts_protein_coding.fa"

#PATH to the vcf file of the ingroup
vcf_file_path="/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/1001g_vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"

#file with ingroup sample IDs that you wish to use from the big vcf file containing 1135 indivudals
pop_info = read.delim('pop_files/iberian_95_pure.txt', header = FALSE)

#importing the conversion table for transitioning from geneIDs - transcriptIDs - peptideIDs
gene_conversion_table = read.delim("input_files/aratha_all_geneIDs_transcriptIDs.txt", header = TRUE, sep = "\t")

#importing then list of geneIDs here
ingroup_geneIDs = read.delim(file = "input_files/aratha_all_geneIDs.txt", header = TRUE)
colnames(ingroup_geneIDs)="geneID"

##### depending on the number of genes, user can define the number of batches and the number of genes in every batch -- ideally a single batch has 10-100 genes ####

## user defined number of batches and number of genes per batch
num_batches = 1 #since the counter starts from 0 - the num_batches will always be total number of batches + 1 
num_genes = 100

#================================================================================================================


#importing relevant packages here
library(ape)
library(stringr)
library(Biostrings)
library(plyr)
library(seqinr)
library(parallel)

#sourcing functions here
source("includes/get_largest_transcript_cds_info.R")
source("includes/get_transcript_reciprocal_alignments.R")
source("includes/get_divergent_sites2.R")
source("includes/get_number_of_synonymous_sites.R")
source("includes/get_nei_gojobori_table.R")
source("includes/get_divergence_stats.R")
source("includes/get_frequencies_cds_for_aratha.R")
source("includes/get_polymorphic_stats.R")
source("includes/polarize_data.R")
source("includes/get_polymorphism_stats_2.R")
source("includes/misc.R")
source("includes/get_syn_nonsyn_positions_poly.R")
source("includes/gff_extraction.R")

#sourcing the main function which contains all the calls
source("main.R")

#declaring constants here
#formality to set threshold of allele frequency within population
maf = 0.0001

#setting min frequency for allele frequency to account for slightly deleterious mutations - workaround for MK test
min_freq = 0.15

#assigning the number of threads for blastn here
num_threads_blastn = 4 #why I use 4 - https://voorloopnul.com/blog/how-to-correctly-speed-up-blast-using-num_threads/

##-- RUN ONLY ONCE -- extracting geneIDs from gff file -- 
##-- depends on the structure of the gff; check corresponding function
#################################################### 

#gff_extraction(ingroup_gff_path)

## done writing file - load this directly ##
####################################################

#importing the gff subset here
ingroup_gff_subset = read.delim("input_files/ingroup_gff_subset.txt", sep = "\t", header = TRUE)
rownames(ingroup_gff_subset) = NULL

#importing the transcript sequences here for both, ingroup and outgroup 
ingroup_transcript_cds = readDNAStringSet(filepath = ingroup_transcript_cds_path, format = "fasta")
outgroup_trasncript_cds = readDNAStringSet(filepath = outgroup_transcript_cds_path, format = "fasta")

#going over all the geneIDs sequentially
for(counter in 0:num_batches){
  
  #counter check
  #counter = 0
  
  #@print check
  print(paste("Going over counter : ",counter,sep = ""))
  
  #declaring start and stop positions
  start_pos = counter*num_genes
  stop_pos = (counter+1)*num_genes
  
  #taking the first batch of geneIDs
  geneID_list = as.character(ingroup_geneIDs[start_pos:stop_pos,])
  
  #creating variable for backup file
  backup_file = paste("backup_",counter,sep="")
  
  #checking if the backup file is created - if not then creating it to store outputs
  if(!dir.exists(paste(backup_file))){
    
    #creating the directory if it does not exist
    dir.create(backup_file)
    
    #print check
    print(paste("Created backup directory : ",backup_file,sep = ""))
    
  }else{print(paste("Found backup directory : ",backup_file,sep=""))}
  
  #giving a call to main.R
  main(geneID_list, backup_file, counter)
  
  #@print check for going over the counter
  print(paste("Gone over counter : ",counter,sep = ""))
}
