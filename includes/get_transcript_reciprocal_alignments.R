#point for making reciprocal blast using local blastn
#the local databases for ingroup and outgroup have already been created 
#these local databases have been created using the cds transcript sequences of both the species

get_transcript_reciprocal_alignments = function(cds_info){

  ## added on 31.10.22: set the blastn DB files as variables, so that they don't contain filename specific arguments
  # previously it was: paste("blastn -db input_files/oligosperma_ref/Coligosperma -query ")
  # now it would be: paste("blastn -db input_files/ingroup_ref/ingroup_prefix ..)
  #ingroup_prefix = unlist(strsplit(list.files('input_files/ingroup_ref', pattern = 'nhr'), '.nhr'))
  #outroup_prefix = unlist(strsplit(list.files('input_files/outgroup_ref', pattern = 'nhr'), '.nhr'))
  
  #declaring transcriptID and seqID
  transcriptID = unique(cds_info$ensembl_transcript_id)
  seqID = as.character(unique(cds_info$seqid))
  geneID = unique(cds_info$ensembl_gene_id)
  
  #@ print check
  print(paste("Extracting sequences, performing reciprocal blast and constructing alignments for :", transcriptID, sep = ""))
  
  ##### forward blastn #####
  
  #importing the transcript sequence from the transcript (cds) sequence file
  transcript_seq = ingroup_transcript_cds[grep(x = names(ingroup_transcript_cds), pattern = transcriptID)]
  
  #storing the transcript sequence to a file for performing lastz alignments
  write.fasta(sequences = transcript_seq, file.out = paste(backup_file,"/transcript_sequences/",
                                                           transcriptID,"_sequence.fa",sep = ""), names = transcriptID)
  
  #forward blast- dummy command
  forward_blastn_main = paste("blastn -db /netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/outgroup_alyra/alyra -query ",backup_file,"/transcript_sequences/",transcriptID,"_sequence.fa ", sep="")
  blastn_fmt = paste("-num_threads ",num_threads_blastn," -outfmt \"6 qseqid qlen qstart qend qseq sseqid slen sstart send sseq bitscore pident evalue\" > ") ##the number of threads assigned to blastn go here
  forward_blastn_name = paste(backup_file,"/transcript_alignments/",transcriptID,"_forward_blast.maf", sep="")
  
  #merging the commands together
  forward_blastn_command = paste(forward_blastn_main, blastn_fmt, forward_blastn_name, sep = "")
  system(forward_blastn_command)
  
  #in case of no alignments, the file is written as null - checking if the written file is not null
  if(file.size(paste(backup_file,"/transcript_alignments/",transcriptID,"_forward_blast.maf", sep="")) > 0){
    
    #reading the forward blast output file
    forward_blastn_output = read.delim(paste(backup_file,"/transcript_alignments/",transcriptID,"_forward_blast.maf", sep=""), header = FALSE, sep = "\t")
    colnames(forward_blastn_output) = c("qseqid", "qlen", "qstart", "qend", "qseq", "sseqid", "slen", "sstart", 
                                        "send", "sseq", "bitscore", "pident", "evalue")
    
    #keeping the column with the highest bitscore only
    forward_blastn_output = forward_blastn_output[forward_blastn_output$bitscore == max(forward_blastn_output$bitscore),]
    
    #keeping only the first entry if there are more than one entries with the highest bitscore 
    ##### ---- keep in mind - this is also the point from where we can detect the duplicated genes - 
    ##### ---- two seperate genes having high bitscore for the same gene
    forward_blastn_output = forward_blastn_output[1,]
    
    #inspired from Urrichio et al 2019 - checking if the identity is atleast 60%
    if(forward_blastn_output$pident > 60){
      
      #storing the geneID and transcriptID of the outgroup top hit
      outgroup_geneID_forward = gsub(pattern = "\\_T.*", replacement = "", x = forward_blastn_output$sseqid)
      outgroup_transcriptID_forward = as.character(forward_blastn_output$sseqid)
      
    }else{
      outgroup_geneID_forward = "NULL1"
      outgroup_transcriptID_forward = "NULL1"
    }
  }else{
    outgroup_geneID_forward = "NULL1"
    outgroup_transcriptID_forward = "NULL1"
  }
  
  ##### reverse blastn #####
  
  #importing the transcript sequence from the transcript (cds) sequence file
  transcript_seq = outgroup_trasncript_cds[grep(x = names(outgroup_trasncript_cds), pattern = outgroup_transcriptID_forward)]
  
  #storing the transcript sequence to a file for performing lastz alignments
  write.fasta(sequences = transcript_seq, file.out = paste(backup_file,"/transcript_sequences/",
                                                           outgroup_transcriptID_forward,"_sequence.fa",sep = ""), names = outgroup_transcriptID_forward)
  
  #reverse command - dummy command -- keeping the fmt command same
  reverse_blastn_main = paste("blastn -db /netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/local_transcript_seq/for_aratha/ingroup_aratha/aratha -query ",backup_file,"/transcript_sequences/",outgroup_transcriptID_forward,"_sequence.fa ", sep="")
  reverse_blastn_name = paste(backup_file,"/transcript_alignments/",transcriptID,"_reverse_blast.maf", sep="")
  
  #merging the commands together
  reverse_blastn_command = paste(reverse_blastn_main, blastn_fmt, reverse_blastn_name, sep = "")
  system(reverse_blastn_command)
  
  #in case of no alignments, the file is written as null - checking if the written file is not null
  if(file.size(paste(backup_file,"/transcript_alignments/",transcriptID,"_reverse_blast.maf", sep="")) > 0){
  
    #reading the reverse blastn output file
    reverse_blastn_output = read.delim(paste(backup_file,"/transcript_alignments/",transcriptID,"_reverse_blast.maf", sep=""), header = FALSE, sep = "\t")
    colnames(reverse_blastn_output) = c("qseqid", "qlen", "qstart", "qend", "qseq", "sseqid", "slen", "sstart", 
                                        "send", "sseq", "bitscore", "pident", "evalue")
    
    #keeping the column with the highest bitscore only
    reverse_blastn_output = reverse_blastn_output[reverse_blastn_output$bitscore == max(reverse_blastn_output$bitscore),]
    
    #keeping only the first entry if there are more than one entries with the highest bitscore 
    ##### ---- keep in mind - this is also the point from where we can detect the duplicated genes - 
    ##### ---- two seperate genes having high bitscore for the same gene
    reverse_blastn_output = reverse_blastn_output[1,]
    
    #inspired from Urrichio et al 2019 - checking if the identity is atleast 60%
    if(reverse_blastn_output$pident > 60){
      
      #extracting the reverse geneIDs
      
      ## this was case specific to the transcript naming in hirsuta, I needed to change this code:
      ## ingroup_geneID_reverse = gsub(pattern = "\\_T.*", replacement = "", reverse_blastn_output$sseqid)
      ## thaliana
      ingroup_geneID_reverse = unlist(strsplit(reverse_blastn_output$sseqid, '[.]'))[1]
      
      ## lyrata
      ## outgroup_geneID_reverse = gsub(pattern = "\\_T.*", replacement = "", reverse_blastn_output$qseqid)
      outgroup_geneID_reverse=reverse_blastn_output$qseqid
      
    }else{
      outgroup_geneID_reverse = "NULL2"
      outgroup_transcriptID_reverse = "NULL2"
      ingroup_geneID_reverse = "NULL3"
    }
  }else{
    outgroup_geneID_reverse = "NULL2"
    outgroup_transcriptID_reverse = "NULL2"
    ingroup_geneID_reverse = "NULL3"
  }
  
  ##### checking condition and processing further #####
  
  # #the ingroup_geneID_reverse has a point - removing this from the geneID - only for humans - this indicates different isoforms
  # ingroup_geneID_reverse = strsplit(ingroup_geneID_reverse, "\\.")[[1]][1]
  
  #putting the bi-condition: a is the best hit for b and b is the best hit for a
  #this condition also makes sure that there is indeed an alignment
  if(outgroup_geneID_forward == outgroup_geneID_reverse & ingroup_geneID_reverse == geneID){
    
    #importing the maf file from the forward blast
    transcript_alignment = read.delim(paste(backup_file,"/transcript_alignments/",transcriptID,"_forward_blast.maf", sep=""), header = FALSE, sep = "\t")
    
    #giving column names to the alignment table
    colnames(transcript_alignment) = c("qseqid", "qlen", "qstart", "qend", "qseq", "sseqid", "slen", "sstart", "send", "sseq", "bitscore", "pident", "evalue")
    
    #selecting the alignment with the highest bitscore
    transcript_alignment = transcript_alignment[transcript_alignment$bitscore == max(transcript_alignment$bitscore),]
    
    #if there are more than one entries with the same bitscore, keeping only the first entry
    transcript_alignment = transcript_alignment[1,]
    
    ### I do not check for percentage identity here as the maf file has already been checked for it during the reciprocal blast step
    
    #checking if there are any N in the sequence alignments, if yes then converting them to "-"
    #subject sequence aka hit sequence
    if(grepl(pattern = "N", ignore.case = TRUE, x = transcript_alignment$sseq)){
      m_seq2 = str_replace_all(string = transcript_alignment$sseq, regex(pattern = "n", ignore_case = TRUE), replacement = "-")
      m_seq2 = unlist(strsplit(m_seq2, split = ""))
    }else{m_seq2 = unlist(strsplit(as.character(transcript_alignment$sseq), split = ""))}
    
    #query sequence
    if(grepl(pattern = "N", ignore.case = TRUE, x = transcript_alignment$qseq)){
      m_seq1 = str_replace_all(string = transcript_alignment$qseq, regex(pattern = "n", ignore_case = TRUE), replacement = "-")
      m_seq1 = unlist(strsplit(m_seq1, split = ""))    
    }else{m_seq1 = unlist(strsplit(as.character(transcript_alignment$qseq), split = ""))}
    
    #making this into a matrix now
    m = matrix(c(m_seq1, m_seq2), nrow = 2, byrow = TRUE)
    
    #checking if there are gaps in the reference matrix
    if("-" %in% m[1,]){
      
      #removing gaps from the reference first
      m = m[,-which((m[1,] == "-"), arr.ind = TRUE)]
      
    }
    
    #making cds positions for the columns names of the matrix
    # listing all of the cds coordinates as the genomic coordinates
    #creating empty list to store it
    cds_pos = c()
    
    #going over all rows sequentially
    for(item in 1:nrow(cds_info)){
      
      #listing a single row here
      row = cds_info[item,]
      
      #all pos in one row
      pos = unlist(row$start):unlist(row$end)
      
      #joining
      cds_pos = c(cds_pos,pos)
      
    }
    
    
    #giving column names to the matrix (on the basis of the length of the aligned region)
    colnames(m) = cds_pos[transcript_alignment$qstart[1]:transcript_alignment$qend[1]]
    
    ###### since the alignments can be out of sync with the open reading frames - I will splice the alignments to be in tune with open reading frame #####
    ## to check if the alignments is in tune with ORF --> check if the start of alignment in query and end of alignment in query is a multiple of 3 -1
    ## so if position is 4 then it is the start of the second codon and if it is 279 then it is the end of 93rd codon
    ##bear in mind the start needs to be multiple of 3 + 1 position and end needs to be a multiple of 3
    if(transcript_alignment$qend %%3 == 0 & transcript_alignment$qstart %% 3 == 1){
      
      #preparing return now
      alignment_matrix = list(alignment = m, transcript_start_pos = transcript_alignment$qstart,
                              id = transcriptID, chr = as.character(unique(cds_info$seqid)))
      
      #return statement
      return(alignment_matrix)
      
    }else{ #if you are inside means you are out of syn with the open reading frame
      
      if(!transcript_alignment$qstart %% 3 == 1){ ## cleaning the start of the transcript
        
        #splicing some positions at the beginning of the alignment
        if(transcript_alignment$qstart %% 3 == 0){
          
          #removing the first column
          m = m[,-1]
          
          #changing the transcript start position
          transcript_alignment$qstart = transcript_alignment$qstart + 1
        }
        
        if(transcript_alignment$qstart %% 3 == 2){
          
          #removing the first two columns
          m = m[,-c(1,2)]
          
          #changing the transcript start position
          transcript_alignment$qstart = transcript_alignment$qstart + 2
          
        }
        
      }
      if(!transcript_alignment$qend %%3 == 0){ ##cleaning the end of the transcript
        
        #splicing some positions at the end of the alignment
        if(transcript_alignment$qend %%3 == 1){
          
          #removing the last column
          m = m[,-ncol(m)]
          
        }
        if(transcript_alignment$qend %%3 == 2){
          
          #removing the last two columns
          m = m[,-c(ncol(m), ncol(m) - 1)]
          
        }
        
      }
      
      #preparing return now
      alignment_matrix = list(alignment = m, transcript_start_pos = transcript_alignment$qstart,
                              id = transcriptID, chr = as.character(unique(cds_info$seqid)))
      
      #return statement
      return(alignment_matrix)
      
    }
    
  }else{return(NULL)}
  
}