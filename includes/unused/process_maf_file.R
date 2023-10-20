#call point for processing the maf file and returning matrix of alignment
process_maf_file = function(path_to_maf_file, cds_info){
  
  #reading lines in the maf file
  maf_file = readLines(path_to_maf_file)
  
  #including only the score lines
  maf_file = maf_file[grepl("s ", maf_file) | grepl("a ", maf_file)]
  
  #removing the first line
  maf_file = maf_file[-1]
  
  #sanity check - checking if there is at all any alignment here
  if(length(maf_file) > 0){
  
    ##block of code for selecting the alignment with the best score and choosing this alignment##
    #going over maf file in blocks of 3 - first entry score, other two are the alignments
  
    #declaring empty dataframe to score the information\
    maf_file_table = data.frame()
    
    #going over list in blocks of 3
    for(item in 0:(length(maf_file)/3 - 1)){
      
      #extracting score line
      score_line = maf_file[(item*3)+1]
      
      #extracting score
      score = str_extract(string = score_line, pattern = "(\\d)+")
      
      #extracting the alignment sequences and its meta-info
      maf_file_subset = gsub(maf_file[(item*3)+(2:3)], pattern = "\\s+", replacement = " ")
      
      #reading the alignment sequences as a dataframe
      df = read.table(text = maf_file_subset, sep = " ")
      
      #removing the first row - it only contains an 's' for sequence
      df = subset(df, select = -c(V1))
      
      #giving column names
      colnames(df) = c("chr", "start", "length", "strand", "stop", "sequence")
      
      #changing datatype
      df$sequence = as.character(df$sequence)
      
      #giving the scores here
      df$score = rep(score)
      
      #binding back to the original dataframe
      maf_file_table = rbind(maf_file_table, df)
    
      }
    
    #changing data-type of score column from character to numeric
    maf_file_table$score = as.numeric(maf_file_table$score)
    
    #subsetting the dataframe to contain the highest scoring alignments only
    maf_file_table = maf_file_table[maf_file_table$score == max(unique(maf_file_table$score)),]
    
    #subsetting the table to contain only the first entry -- just in case two or more alignments
    #contain the same high score
    maf_file_table = maf_file_table[1:2,]
    
    #getting sequence information - first entry is always the outgroup and second in the species of interest
    m_seq1 = unlist(strsplit(maf_file_table$sequence[2], split = ""))
    m_seq2 = unlist(strsplit(maf_file_table$sequence[1], split = ""))
    
    ## converting cds coordinates to their genomic coordinates ##
    # listing all of the cds coordinates as the genomic coordinates
    #creating empty list to store it
    cds_pos = c()
    
    #going over all rows sequentially
    for(item in 1:nrow(cds_info)){
      
      #listing a single row here
      row = cds_info[item,]
      
      #all pos in one row
      pos = row$start:row$end
      
      #joining
      cds_pos = c(cds_pos,pos)
      
    }
    
    #converting into a matrix
    m = matrix(c(m_seq1, m_seq2), nrow = 2, byrow = TRUE)
    
    #checking if there are gaps in the reference matrix
    if("-" %in% m[1,]){
    
      #removing gaps from the reference first
      m = m[,-which((m[1,] == "-"), arr.ind = TRUE)]
    
    }
    
    #giving the column names now
    colnames(m) = cds_pos[1:dim(m)[2]]
    
    ## adding the return of the transcript start position to facilitate syn and nonsyn prediction later
    transcript_start_pos = maf_file_table$start[2]
    
    #making list to return
    
    return(list(m, transcript_start_pos))
  }else{return(NULL)}
  
  }