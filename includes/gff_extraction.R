## gff file manipulation: 
## extracting geneIDs 
## preparing a subset file that would contain only genes (that have cds information) and cds 

## NOTE:
##-- depends on the structure of the organism's gff and naming patterns
## for Thaliana:
##               genes = AT1G17000
##               cds=    AT1G17000.1
## and inside the gff the attributes colum has the following structure
## "ID=gene:AT1G17000;Name=ATTPS3;biotype=nontranslating_CDS;description=trehalose-phosphatase/synthase 3 [Source:TAIR%3BAcc:AT1G17000];gene_id=AT1G17000;logic_name=araport11"

## for Humans: instead of having a geneID and transcriptID; humans have geneID, peptideID and transcriptID
##              genes         =       ENSG00..186092
##              peptides      =       ENSP00...2345
##              transcripts   =       ENST00....1234

gff_extraction=function(ingroup_gff_path) {

#reading the gff file
ingroup_gff = read.gff(ingroup_gff_path, GFF3 = TRUE)

#--- exploratory part ----------
# find genes that do not have cds information

#keep gene IDS:
genes = ingroup_gff[ingroup_gff$type %in% c("gene"),"attributes"]
#-- the atributes field is a big string, we are only keeping the ID field:
genes = gsub(".*ID=gene:", "" , genes) 
#-- to remove the tail of the strings (contains extra fields)
genes = sub(';.*', '',genes)

# # keep cds IDS: --- this does not apply for humans as humans do not just have a geneID and transcriptID but a geneID, transcriptID and peptideID
# cds= ingroup_gff[ingroup_gff$type %in% c("CDS"),"attributes"]
# #-- the atributes field is a big string, we are only keeping the ID field:
# cds = gsub(".*ID=CDS:", "" , cds) 
# #-- to remove the tail of the strings (contains extra fields)
# cds = sub(';.*', '',cds)

## exploring the peptide IDs and transcript IDs - no need to run it this is just for exploration
#accessing the CDS information first
cds= ingroup_gff[ingroup_gff$type %in% c("CDS"),"attributes"]

#for peptides
#-- the atributes field is a big string, we are only keeping the ID field:
peptides = gsub(".*ID=CDS:", "" , cds)
#-- to remove the tail of the strings (contains extra fields)
peptides = sub(';.*', '',peptides)

#for transcripts
#-- the atributes field is a big string, we are only keeping the ID field:
transcripts = gsub(".*Parent=transcript:", "" , cds)
#-- to remove the tail of the strings (contains extra fields)
transcripts = sub(';.*', '',transcripts)

# # From each cds ID, extract what is the gene name that it corresponts to. 
# # Keep unique (because >1 cds can correspont to the same gene ID)
# genes_corresponding_to_cds=unique(unlist(lapply(strsplit(cds, '[.]'), '[[', 1)))
# 
# #SO: 28 genes DO NOT have CDS information in the gff
# length(genes)
# length(which(genes %in% genes_corresponding_to_cds))
# 
# # those genes are: 
# orphan_genes= setdiff(genes, genes_corresponding_to_cds)

#---- Those genes would be removed for any downsteam analysis ------


#subsetting the gff file to only contain genes and CDS
ingroup_gff = ingroup_gff[ingroup_gff$type %in% c("gene", "CDS"),]

#-- the atributes field is a big string, we are only keeping the ID field:
#-- if it is a gene:
ingroup_gff$attributes = gsub(".*ID=gene:", "" , ingroup_gff$attributes) 
#-- if it is a CDS:
ingroup_gff$attributes = gsub(".*ID=CDS:", "" , ingroup_gff$attributes) 

#-- to remove the tail of the strings (contains extra fields)
ingroup_gff$attributes= sub(';.*', '',ingroup_gff$attributes)

#### for humans - the choice of geneIDs is already made and the conversion table of the selected transcript (peptides is stored in input_file)
#### for more read -- input_files/readme

# #listing the geneIDs here
# ingroup_gff_geneIDs = ingroup_gff$attributes[ingroup_gff$type == "gene"]
# 
# #-- removing genes with no cds information aka 'the orphan genes'
# ingroup_gff_geneIDs=setdiff(ingroup_gff_geneIDs, orphan_genes)
# 
# 
# #writing list of geneIDs to a file
# ingroup_gff_geneIDs = data.frame(geneID = ingroup_gff_geneIDs)
# write.table(file = "input_files/ingroup_geneIDs.txt", x = ingroup_gff_geneIDs$geneID, col.names = TRUE,
#             row.names = FALSE, quote = FALSE)

#writing the subset of gff to a file
write.table(file = "input_files/ingroup_gff_subset.txt", x = ingroup_gff,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


}