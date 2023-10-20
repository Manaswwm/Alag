convert_cds_to_genomic_coordinates=function(x){
  

  server <- "http://rest.ensembl.org"
  ext=paste("/map/cds/",x[1],"/",gsub(" ","",x[2]),"..",gsub(" ","",x[2]),"?include_original_region=1", sep="")
  
  print(ext)
#  print(paste(x[1],x[2]))

  r <- GET(paste(server, ext, sep = ""), content_type("application/json"), config = httr_config)
  
  response=fromJSON(toJSON(content(r), force = TRUE))
#  print(is.null(names(response$mappings$original)))
  
  if(r$status_code==200 && dim(response$mappings$original)[2]!=0){#BIG UNCERTAINTY HERE, WHY sSECOND condition-> mention to ensembl helpdesk
    
#  print(names(response$mappings$original))
  names(response$mappings$original)<-paste("original",names(response$mappings$original), sep=".")
  names(response$mappings$mapped)<-paste("mapped",names(response$mappings$mapped), sep=".")

  response=cbind(response$mappings$original, response$mappings$mapped)
  
  
  return(response)
  
  } else return(NULL) }