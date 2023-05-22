suppressWarnings ( suppressPackageStartupMessages ( library ( data.table ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( dplyr ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( fst ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( ggplot2 )))
library(dndscv)


#1
###################################################################################
################## Somatic mutation data 【mutation table】########################

mut_table = read_fst('figurecode/data/20221109_TRACERx421_mutation_table.fst')
mutation <- mut_table[,c("tumour_id","chr","start","ref","var")]
# rename columns to match ones in dndscv package
colnames(mutation) = gsub('start', 'pos', colnames(mutation))

###### !!! to make chr readable for dndscvpackage(chr column only number)
mutation$chr = gsub("chr","",as.vector(mutation$chr))
head(mutation)
# tumour_id  chr      pos ref var
#1  CRUK0005 19 34291428   G   T
#2  CRUK0005 19  9084770   C   A
#3  CRUK0005  1  2160418   G   T
rm(mut_table)
################################################################################





#2
################################################################################
################ Weights for mutational signatures ##########################

#(1)
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
signature_weights_clonalMuts_perTumour <- signature_weights[["signature_weights_clonalMuts_perTumour"]]
signature_weights_subclonalMuts_perTumour <- signature_weights[["signature_weights_subclonalMuts_perTumour"]]
signature_weights_privateMuts_perTumour <- signature_weights[["signature_weights_privateMuts_perTumour"]]

#(2)
############ signature weights group tumour_id ################

number <- c(1, 2, 4, 5, 13, 92, 44)

#define signature weights<0.1 as undetected
SBS_list <- list()
for (i in number) {
  col <- paste0("SBS", i)
  filter_table <- signature_weights_perTumour[[col]] > 0.1
  table <- subset(signature_weights_perTumour, filter_table, select = c("tumour_id"))
  SBS_list[[col]] <- table
}

for (i in 1:length(SBS_list)) {
  table <- SBS_list[[i]]
  merged_table <- merge(mutation,table,  by = "tumour_id")
  SBS_list[[i]] <- merged_table
}

###############################################################################




#2
###############################################################################
###################### run dndscv on SBS mutations ############################

name <- c("1", "2", "4", "5", "13", "92", "44")

dndsout_SBS <- list()
for (i in name){
  SBSi <- paste0("SBS", i)
  sublist <- SBS_list [[SBSi]]
  dndsout <- dndscv(sublist, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)
  dndsout_SBS[[i]] <- assign(paste0("dndsout_", SBSi), dndsout, envir = .GlobalEnv)
}  


#(1) 
############# output significant genes ##############
signature_driver_q0.1_list <- list()
for (i in name) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBSi <- get(table)
  sel_cv_subtype <- dndsout_SBSi$sel_cv
  signature_driver_q0.1_list[[i]] <- sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")]
}

for (i in name) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBSi <- get(table)
  sel_cv_subtype <- dndsout_SBSi$sel_cv
  print(head(sel_cv_subtype), digits = 3) 
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}  

for (i in name) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv","wmis_cv","wnon_cv","wspl_cv","wind_cv")], digits = 3)
}  


#(2) 
########## output global dN/dS Global dN/dS estimates [globaldnds] ############

for (i in name) {
  print(i)
  print("Global dN/dS")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)
  print (dndsout_SBS$globaldnds)
}  


#(3)
######### theta/θ value to test suitability for dNdScv #######

for (i in name) {
  print(i)
  print("theta value")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  print(dndsout_SBS$nbreg$theta)
}


#(4)
############### AIC model to measure the fit of a statistical model #############
# the smaller the AIC, the better the model

for (i in name) {
  print(i)
  print("AIC")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  print(AIC(dndsout_SBS$poissmodel))
}

########################################################################
# 针对warning:1: In dndscv(sublist, max_muts_per_gene_per_sample = Inf,  ... :
# Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.
# To check with this sites:
mutations$pos <- as.numeric(mutations$pos)
mutations <- mutations[order(mutations$tumour_id, mutations$chr, mutations$pos),]
ind = which(diff(mutations$pos)==1)
mutations[unique(sort(c(ind,ind+1))),]
