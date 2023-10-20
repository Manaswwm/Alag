#call point for getting the transcript alignments given cds coordinates
get_transcript_alignments = function(cds_info){
  
  #declaring transcriptID and seqID
  transcriptID = unique(cds_info$transcriptID)
  seqID = as.character(unique(cds_info$seqid))
  
  #@ print check
  print(paste("Extracting sequences for :", transcriptID, sep = ""))
  
  #importing the transcript sequence from the transcript (cds) sequence file given by Christos
  transcript_seq = hirsuta_transcript_cds[names(hirsuta_transcript_cds) == transcriptID]
  
  #storing the transcript sequence to a file for performing lastz alignments
  write.fasta(sequences = transcript_seq, file.out = paste(backup_file,"/transcript_sequences/",
                            transcriptID,"_sequence.fa",sep = ""), names = transcriptID)
  
  ## LASTZ block ##
  
  #@ print check
  print(paste("Making alignments for :", transcriptID, sep = ""))
  
  #importing the sequnece files  
  lastz_seq_info = paste("lastz input_files/oligosperma_reference/C.oligosperma.fa[multiple] ",
                         backup_file,"/transcript_sequences/",transcriptID,"_sequence.fa ", sep = "")
  
  #I allow for ambiguous characters in the sequence file, I get only one hit with HSPBEST = 1
  #I get the output in maf format so that it is easier for me to process sequences later
  lastz_misc_info = paste("--ambiguous=iupac --queryhspbest=1 --format=maf --strand=plus ")

  #declaring the output file name
  lastz_output = paste("--output=",backup_file,"/transcript_alignment/",transcriptID,"_alignment.maf",sep="")
  
  #putting things together here
  lastz_command = paste(lastz_seq_info, lastz_misc_info, lastz_output, sep="")
  
  #print check
  #print(paste("Performing alignment of entire transcript sequence (merged) of : ",transcriptID,sep=""))
  
  #calling lastz here and storing alignment directly
  system(lastz_command)
  
  #print check
  print(paste("Done making alignment of transcript ID : ",transcriptID, sep = ""))
  
  ##LASTZ block over - processing the maf files now##
  
  #@ print check
  print(paste("Processing alignments for :", transcriptID, sep = ""))
  
  #giving the path to the written maf file to extract the substitution matrix
  alignment_matrix = process_maf_file(paste(backup_file,"/transcript_alignment/",transcriptID,"_alignment.maf",sep=""),
                                      cds_info)
  
  #adding back the transcriptID and chr to the matrix
  alignment_matrix=list(alignment=alignment_matrix[[1]], transcript_start_pos = alignment_matrix[[2]],
                        id=transcriptID, chr=seqID)
  
  #returninng the alignment
  return(alignment_matrix)
  
  
}

# ### dummy input -- remove before launching for a batch!! ###
# #cds_info_sub = cds_info[[1]]
# 
# #getting the sequence of the entire transcript from cds coordinates
# #declaring an empty string to store the transcript sequences
# transcript_seq = ""
# 
# #declaring transcriptID
# transcriptID = as.character(unique(cds_info$transcriptID))
# 
# #print check
# print(paste("Making transcript alignments for id : ", transcriptID, sep = ""))
# 
# #declaring seqid
# seqID = unique(as.character(cds_info$seqid))
# 
# #declaring the fasta sequence file from which I want to extract the sequence
# fa_file = readDNAStringSet(paste("dummy_input_files/Athaliana_seq/Athaliana_",seqID,".fa", sep =""))
# 
# #going over all start and end positions
# for(item in 1:nrow(cds_info)){
#   
#   #one row in a dataframe
#   row = cds_info[item,]
#   
#   #subsetting the chromosome fasta to get the cds_sequence
#   cds_seq = toString(subseq(fa_file, start = row$start, end = row$end))
#   
#   #merging info from all the rows here
#   transcript_seq = paste(transcript_seq, cds_seq, sep = "")
#   
# }
# 
# #checking if this transcript is on a negative strand - if yes then reverse complementing this
# #sequence before the output
# if(unique(cds_info$strand) == "-"){
#   
#   #first converting this in to a DNA string
#   transcript_seq = DNAString(x = transcript_seq)  
#   
#   #now reverse complementing
#   transcript_seq = as.character(reverseComplement(transcript_seq))
#   
# }
# 
# 
# #storing the transcript sequence to a file for performing lastz alignments
# write.fasta(sequences =  transcript_seq, file.out = paste(backup_file,"/transcript_sequences/",
#                           transcriptID,"_sequence.fa",sep = ""), names = transcriptID)
# 