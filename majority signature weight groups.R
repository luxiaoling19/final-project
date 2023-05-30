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
################ highest Weights for mutational signatures ##########################

#(1)
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
rm(signature_weights)

#(2) only include the highest weights
############ signature weights group tumour_id ################
tumour_id <- signature_weights_perTumour[["tumour_id"]]
max_SBS <- data.table(tumour_id = character(), max_sig_weights = character())  # create result table

for (i in 1:nrow(signature_weights_perTumour)) {
  row <- signature_weights_perTumour[i, ]
  max_value <- max(na.omit(as.numeric(row[-1])))    # ignore column 1(tumour_id)
  column_index <- which(row[-1] == max_value)
  column_name <- colnames(signature_weights_perTumour)[column_index + 1] # the first column is ignored, the index needs to be added by 1
  max_SBS <- rbind(max_SBS, data.table(tumour_id = tumour_id[i], max_sig_weights = column_name))
}

head(max_SBS)

###############################################################################




#3
######################################################################################
################  mutational signatures group by maximum weights ###################
subtype <- c("1","2", "4", "5", "13", "92", "44")

max_weights_sig_list <- list()
for (i in subtype){
  SBSi <- paste0("SBS", i)
  max_weights_tumour_id <- max_SBS[max_sig_weights==SBSi][,c("tumour_id")]
  max_weights_sig_list[[SBSi]] <- max_weights_tumour_id
}


for (i in 1:length(subtype)) {
  SBSi <- paste0("SBS", subtype[i])
  table <- max_weights_sig_list[[SBSi]]
  merged_table <- merge(mutation, table, by = "tumour_id")
  max_weights_sig_list[[SBSi]] <- merged_table
}


head(max_weights_sig_list)
######################################################################################




#4
###############################################################################
###################### run dndscv on SBS mutations ############################
subtype <- c("2", "4", "5", "13", "92", "44")
dndsout_max_SBS <- list()
for (i in subtype){
  SBSi <- paste0("SBS", i)
  sublist <- max_weights_sig_list[[SBSi]]
  dndsout <- dndscv(sublist, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)
  dndsout_max_SBS[[SBSi]] <- assign(paste0("dndsout_", SBSi), dndsout, envir = .GlobalEnv)
}  


#(1) 
############# output significant genes ##############

for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBSi <- get(table)
  sel_cv_SBSi <- dndsout_SBSi$sel_cv
  print(sel_cv_SBSi[sel_cv_SBSi$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")])
}

for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBSi <- get(table)
  sel_cv_SBSi <- dndsout_SBSi$sel_cv
  print(head(sel_cv_SBSi), digits = 3) 
  print(sel_cv_subtype[sel_cv_SBSi$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}  



#(2) 
########## output global dN/dS Global dN/dS estimates [globaldnds] ############

for (i in subtype) {
  print(i)
  print("Global dN/dS")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)
  print (dndsout_SBS$globaldnds)
}  


#(3)
######### theta/θ value to test suitability for dNdScv #######

for (i in subtype) {
  print(i)
  print("theta value")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  print(dndsout_SBS$nbreg$theta)
}


#(4)
############### AIC model to measure the fit of a statistical model #############
# the smaller the AIC, the better the model

for (i in subtype) {
  print(i)
  print("AIC")
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  print(AIC(dndsout_SBS$poissmodel))
}


######## other outputs ########
head(dndsout_ASC$annotmuts)
# 【annotmuts】an annotated table of coding mutations编码突变的注释表 
# 【mle_submodel】MLEs of mutation rate parameters突变率参数的最大似然估计 
# 【genemuts】 a table with the observed expected number of mutations per gene 每个基因观察到的和预期的突变数量表
annotmuts_major_weights <- list()
for (i in subtype) {
  print(i)
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  annotmuts_major_weights[[table]] <- dndsout_SBS$annotmuts
}
genemuts_major_weights <- list()
for (i in subtype) {
  print(i)
  table <- paste0("dndsout_SBS", i) 
  dndsout_SBS <- get(table)  
  genemuts_major_weights[[table]] <- dndsout_SBS$genemuts
}

library(writexl)
write_xlsx(annotmuts_major_weights, "annotmuts_major_weights.xlsx")
write_xlsx(genemuts_major_weights, "genemuts_major_weights.xlsx")



