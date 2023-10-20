#short script that takes imports the transcript sequences file per species and subsets to keep only protein
#coding transcript sequences

##### setting paths for the input and output transcript sequences to the ingroup and the outgroup species ######

#ingroup
input_ingroup_transcript_sequences = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/ingroup_aratha/Arabidopsis_thaliana.TAIR10.cds.all.fa"
output_ingroup_transcript_sequences = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/ingroup_aratha/aratha_transcripts_proteincoding.fa"


#outgroup
input_outgroup_transcript_sequences = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/outgroup_alyra/Arabidopsis_lyrata.v.1.0.cds.all.fa"
output_outgroup_transcript_sequences = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/outgroup_alyra/alyra_transcripts_protein_coding.fa"

##### done setting paths #####

#importing relevant packages
library(Biostrings)

#first for thaliana
aratha_transcript_cds = readDNAStringSet(filepath = input_ingroup_transcript_sequences, format = "fasta")

#subsetting to keep only protein coding transcripts
aratha_transcript_cds_proteincoding = aratha_transcript_cds[grep(x = names(aratha_transcript_cds), pattern = "protein_coding")]

#writing back to the same file
writeXStringSet(x = aratha_transcript_cds_proteincoding, filepath = output_ingroup_transcript_sequences, format = "fasta")

#second for lyrata
alyra_transcrip_cds = readDNAStringSet(filepath = input_outgroup_transcript_sequences, format = "fasta")

#keeping only protein coding transcripts
alyra_transcript_cds_proteincoding = alyra_transcrip_cds[grep(x = names(alyra_transcrip_cds), pattern = "protein_coding")]

#writing back to the same folder
writeXStringSet(x = alyra_transcript_cds_proteincoding, filepath = output_outgroup_transcript_sequences, format ="fasta")
