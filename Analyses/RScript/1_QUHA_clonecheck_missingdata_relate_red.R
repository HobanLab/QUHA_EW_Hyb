##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

##############################################
############# Convert to Genind ##############
##############################################
setwd("C:\\Users\\eschumacher\\Documents\\Qhavardii_ex_situ-main")

QUHA_wild_gen <- read.genepop("QH_wild.gen", ncode = 3)


################################################################ 
########## Run Clone Check and Relatedness Analysis ############
################################################################
####run clone check 
##convert to genelcone object
#QUAC_geneclone <- as.genclone(QUAC_garden_wild_gen)

QUHA_geneclone <- as.genclone(QUHA_wild_gen)

##identify multi-locus genotypes (non-clones)
#QUAC_gen_id <- mlg.id(QUAC_geneclone)

QUHA_mlg_id <- mlg.id(QUHA_geneclone)

##function to pull out all clones into a list
QUHA_clone_index <- which(sapply(QUHA_mlg_id, function(x) length(x)>1))

##save a data frame of clones
QUHA_clone_list <- list()

##write a loop to compare all clones 
for(i in 1:length(QUHA_clone_index)){
  
  QUHA_clone_list[[i]] <- QUHA_mlg_id[[paste0(QUHA_clone_index[[i]])]]
  
}

##now remove clones from the matrix 
QUHA_noclones <- clonecorrect(QUHA_geneclone)
##convert back to a genind object
QUHA_genind_nocl <- genclone2genind(QUHA_noclones) ##left with 542 individuals
##remove individuals with too much missing data 
QUHA_genind_nomd <- missingno(QUHA_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) ##there were 7 individuals with too much missing data, left with 449 individuals

##write out to a genalex
genind2genalex(QUHA_genind_nomd, "C:\\Users\\eschumacher\\Documents\\GitHub\\QUHA_EW_Hyb\\QUHA_noclone_nomd_genalex.csv")

##read in relatedness file 
QUHA_relate_df <- read.csv("C:\\Users\\eschumacher\\Documents\\GitHub\\QUHA_EW_Hyb\\Data_Files\\QUHA_noclone_nomd_df.csv")

##run relatedness analysis on the cleaned score data frame 
QUHA_relatedness_df <- Demerelate(QUHA_relate_df, object = T, value = "loiselle")
##now identify how many individuals have greater than 25% relatedness = half siblings
QUHA_halfsib_names <- names(which(unlist(QUHA_relatedness_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
#QUAC_halfsib_names_cleanfront <- gsub("^.*\\.","", QUAC_halfsib_names)
QUHA_halfsib_names_cleanfront <- gsub("^.*\\.","", QUHA_halfsib_names)

#QUAC_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_halfsib_names_cleanfront)
QUHA_halfsib_names_cleanback <- gsub("^.*\\_","", QUHA_halfsib_names_cleanfront)

#QUAC_relate_ind_remove <- unique(QUAC_halfsib_names_cleanback) ##260 of these individuals have > 25% relatedness
QUHA_relate_ind_remove <- unique(QUHA_halfsib_names_cleanback) ##252 individuals are 25% related and > 

##now limit genind object by these names 
#QUAC_relate_red_gen <- QUAC_genind_nomd[!rownames(QUAC_genind_nomd@tab) %in% QUAC_relate_ind_remove,]
##and export to genalex
#genind2genalex(QUAC_relate_red_gen, file="QUAC_data_frames/Relate_Red/QUAC_relate_red_garden_wild_genalex.csv",
#overwrite = TRUE)
##limit data frame with the indiviudals that are highly related  
#QUAC_red_relate_garden_wild_df <- QUAC_rel_df[!QUAC_rel_df$ID %in% QUAC_relate_ind_remove,]



