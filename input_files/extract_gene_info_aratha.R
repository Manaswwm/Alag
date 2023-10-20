#short code to obtain all the possible geneIDs within aratha - for Chr1:Chr5

#importing relevant packages
library(biomaRt)

#declaring mart for thaliana
thaliana_mart = useMart(host="https://plants.ensembl.org", "plants_mart",
                        dataset = "athaliana_eg_gene")

#getting geneIDs and corresponding information on the transcripts per gene
geneIDs_aratha = getBM(attributes = c("ensembl_gene_id", "chromosome_name", "gene_biotype",
                                    "ensembl_transcript_id", "transcript_is_canonical"), filters = "chromosome_name",
                     values = seq(1, 5, by = 1), useCache = FALSE, mart = thaliana_mart)

#first subsetting the dataframe to contain only protein_coding genes
geneIDs_aratha_subset = geneIDs_aratha[geneIDs_aratha$gene_biotype == "protein_coding",]

#### unlike hsap where we choose MANE transcript - for aratha we will choose canonical transcript ####
geneIDs_aratha_subset = geneIDs_aratha_subset[!is.na(geneIDs_aratha_subset$transcript_is_canonical),]

#double checking that the geneIDs are not repeated
length(unique(geneIDs_aratha_subset$ensembl_gene_id)) #not repeated - so every row is a gene and its corresponding "chosen" transcript

#saving only the geneIDs column first
write.table(file = "aratha_all_geneIDs.txt", x = geneIDs_aratha_subset$ensembl_gene_id, col.names = TRUE, row.names = FALSE,
            quote = FALSE)

#next saving the transcript ID and the geneID for choosing the transcript in alag
write.table(file = "aratha_all_geneIDs_transcriptIDs.txt", x = geneIDs_aratha_subset[,c("ensembl_gene_id", "ensembl_transcript_id")],
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
