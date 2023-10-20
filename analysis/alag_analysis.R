##the script reads output files from Alag, processes them and creates outputs (constraint ratios, alpha estimates, SFS visualization)

#importing relevant libraries
library(ggplot2)
library(dplyr)

######################################################################
# ------------------- setting paths here ------------------- #

#path to the folder containing the output files named --> poly_div_gene_*
poly_div_stats_folder_path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_IB/"

#path to the folder containing the backup folders per batch --> backup_*
backup_files_path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_IB/"

######################################################################

#### Step 0 - reading the input files ####

#listing filenames containing the poly-div stats per batch
poly_div_stats_files = list.files(path = poly_div_stats_folder_path, pattern = "poly_div_gene_*", full.names = TRUE)

#aggregating the listed files to form a single table 
poly_div_stats = lapply(poly_div_stats_files, function(x){read.delim(x, sep = "\t")})
poly_div_stats = do.call(rbind, poly_div_stats)
poly_div_stats = unique(poly_div_stats) ##in some cases batch IDs overlap - removing them here

#listing backup folders per batch
backup_files = list.files(path = backup_files_path, pattern = "backup_*", full.names = TRUE)

#declaring an empty dataframe to fill in the information on frequencies 
freq_cds = data.frame()

#going over all batch files sequentially
for(file in backup_files){
  
  #reading a single batch file
  load(paste(file,"/frequencies_cds", sep = ""))
  
  #merging
  freq_cds = rbind(freq_cds, frequencies_cds)
  
  
}

#taking unique
freq_cds = unique(freq_cds) #ensuring no point mutations are repeating due to batch overlapping genes

#cleaning up
rm(poly_div_stats_files, poly_div_stats_folder_path, backup_files_path, file, backup_files, frequencies_cds)


#### Step 1 - Calculating the constraint ratios ####

#estimating polymorphism statistics
pi_nonsyn = mean(poly_div_stats$pi_nonsyn)
pi_nonsense = mean(poly_div_stats$pi_nonsense)
pi_syn = mean(poly_div_stats$pi_syn)
pi_n_pi_s = pi_nonsyn/pi_syn
pi_nonsense_pi_s = pi_nonsense/pi_syn

#estimating divergence statistics
div_nonsyn = mean(poly_div_stats$div_nonsyn)
div_nonsense = mean(poly_div_stats$div_nonsense)
div_syn = mean(poly_div_stats$div_syn)
k_n_k_s = div_nonsyn/div_syn
k_nonsense_k_s = div_nonsense/div_syn

#aggregating polymorphism and divergence statistics in a dataframe
constraint_ratios = data.frame(pi_syn = pi_syn, pi_nonsyn = pi_nonsyn, pi_n_pi_s = pi_n_pi_s, pi_nonsense = pi_nonsense, div_syn = div_syn, 
                               div_nonsyn = div_nonsyn, k_n_k_s = k_n_k_s, div_nonsense = div_nonsense, k_nonsense_k_s = k_nonsense_k_s)

## writing the polymorphism and divergence statistics to an output file
write.table(file = "constraint_ratios.txt", x = constraint_ratios, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#cleaning up
rm(constraint_ratios, pi_nonsyn, pi_nonsense, pi_syn, pi_n_pi_s, pi_nonsense_pi_s, div_syn, div_nonsyn, div_nonsense, k_n_k_s, k_nonsense_k_s)

#### Step 2 - Constructing site frequency spectrum for nonsyn and syn variants occurring within the population #####

#creating a seperate dataframe to contain the nonsyn and syn frequencies
freq_info = freq_cds[,c("freq_der", "type")]

#discretizing the frequencies in 10 classes
freq_info$class = cut(freq_info$freq_der, breaks = seq(0, 1, by=0.1))

#constructing sfs using dplyr function
sfs = freq_info %>%
  rowwise() %>%
  group_by(type, class) %>%
  tally(name = "n") %>%
  mutate(prop = n/sum(n))
  
#adding another column for frequency to be plotted in the sfs
sfs$freq = rep(seq(0.1,1,by = 0.1),3)

#plotting the sfs
ggplot(sfs, aes(freq, y = prop, fill = type)) + geom_bar(stat = "identity", position = position_dodge(), colour = "black", alpha = 0.8) + theme_classic()+
  xlab("Frequencies") + ylab("Proportion") + labs(fill = "Type") + scale_x_continuous(breaks = seq(0,1,by = 0.1))

#saving the output plot to a file
ggsave(filename = "sfs_aratha_iberia.png", dpi= 300, device = "png")

#cleanin up
rm(sfs, freq_info)

#### Step 3 - Constructing class-specific alpha estimates ####

#splitting the frequencies information for nonsyn and syn
freq_info_nonsyn = freq_cds[freq_cds$type == "nonsyn",]
freq_info_syn = freq_cds[freq_cds$type == "syn",]

## ! declaring mean values on lengths and the divergent time scale that will be used the same for every class in the SFS ! ##
div_syn = mean(poly_div_stats$div_syn)
div_nonsyn = mean(poly_div_stats$div_nonsyn)
nonsyn_length = sum(poly_div_stats$nonsyn_length)
syn_length = sum(poly_div_stats$syn_length)

#### !!! done declaring constants !! ####

#defining an empty dataframe to fill in the alpha values per sfs class
alpha_sfs = data.frame()

#going over all sfs classes sequentially - per class is 0.01
for(int in seq(0.01, 0.99 ,by = 0.01)){
  
  #extracting the frequencies belonging to these intervals
  syn_freq = freq_info_syn$freq_der[freq_info_syn$freq_der >= int & freq_info_syn$freq_der < (int+0.01)]
  nonsyn_freq = freq_info_nonsyn$freq_der[freq_info_nonsyn$freq_der >= int & freq_info_nonsyn$freq_der < (int+0.01)]
  
  #calculating pi
  pi_s = sum(2 * syn_freq * (1 - syn_freq))/syn_length
  pi_n = sum(2* nonsyn_freq * (1 - nonsyn_freq))/nonsyn_length
  
  #calculating alpha
  alpha = 1- ((div_syn/div_nonsyn) * (pi_n/pi_s))
  
  #binding together
  df = data.frame(alpha = alpha, int = int)
  
  #merging
  alpha_sfs = rbind(alpha_sfs, df)
  
}

#removing NAN estimates
alpha_sfs = alpha_sfs[!is.nan(alpha_sfs$alpha),]

#plotting the distribution of alpha per sfs class
ggplot(data = alpha_sfs, aes(x = int, y = alpha)) + geom_point(size = 4) + scale_x_continuous(breaks = seq(0,1, by = 0.1)) +theme_bw() + 
  xlab("SFS class") + ylab(bquote(~"\u03B1 estimates")) + ggtitle(bquote("Distribution of"~"\u03B1 estimates per frequency class"))+
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)) + scale_y_continuous(breaks = seq(-0.9, 0.4, by = 0.1)) #6 rows are removed due to the low pop size - read comments on lines 766 and 767

#saving the plot to a file
ggsave(filename = "wgs_alpha_estimates_aratha_IB.png", device = "png", width = 10, height = 10)