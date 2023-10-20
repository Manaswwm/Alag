# call point for annotating polymorphic variants from transcript specific vcf file
annotate_polymorphic_variants = function(path_to_vcf_folder, transcriptID){
  
  #checking if there is a vcf file corresponding to the transcript ID
  if(file.exists(paste(path_to_vcf_folder,transcriptID,".vcf",sep=""))){
    
    #print check
    print(paste("Annotating variants for transcript ID : ",transcriptID,sep = ""))
    
    #constructing parts of the snpEff command
    #main command body
    snpeff_main = paste("java -Xmx4g -jar ../../tools/snpEff/snpEff.jar CarHir_2 ")
    
    #input and output vcf files
    snpeff_vcf_files = paste(path_to_vcf_folder,transcriptID,".vcf > ",path_to_vcf_folder,transcriptID,"_annotated.vcf", sep="")
    
    #print check
    print(paste(snpeff_main,snpeff_vcf_files,sep = ""))
    
    #submitting to system
    system(paste(snpeff_main,snpeff_vcf_files,sep=""),intern = TRUE)
    
  } #done annotating variants for this transcript
}