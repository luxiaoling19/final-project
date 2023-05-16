suppressWarnings ( suppressPackageStartupMessages ( library ( data.table ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( dplyr ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( fst ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( ggplot2 )))


#1
###################################################################################
################## Somatic mutation data 【mutation table】########################

mut_table = read_fst('figurecode/data/20221109_TRACERx421_mutation_table.fst')
mutation <- mut_table[,c("tumour_id","chr","start","ref","var")]
# rename columns to match ones in dndscv package
colnames(mutation) = gsub('start', 'pos', colnames(mutation))

###### !!! to make column "chr" readable for dndscvpackage(chr column only number)
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
################# Pathology data 【Tumour Histological subtype】 ######################

tumor_df = readRDS('figurecode/data/20221109_TRACERx421_all_tumour_df.rds')
# rename columns to match ones in mutation table and CNA
colnames(tumor_df) = gsub('tumour_id_muttable_cruk', 'tumour_id', colnames(tumor_df))
colnames(tumor_df) = gsub('cruk_id', 'patient_id', colnames(tumor_df))
head(tumor_df)
########## Overview plots ##########
## patients' distribution across tumor subtypes(Histological subtype)
patient_by_subtype = as.data.table(tumor_df)
patient_by_subtype = patient_by_subtype[,(length(unique(tumour_id))),
                                        by = Histology_per_tumour_id_muttable]
colnames(patient_by_subtype) = c('subtype', 'n')
patient_by_subtype = patient_by_subtype[order(patient_by_subtype$n)]
patient_by_subtype$subtype = factor(patient_by_subtype$subtype,
                                    patient_by_subtype$subtype)
ggplot(data = patient_by_subtype, 
       aes(x = subtype, y = n, fill = subtype, label = n)) +
  geom_bar(stat = "identity") + 
  geom_text(vjust = -1, size = 3) +
  xlab('Histological subtype') + ylab('Number of tumors') +
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
       legend.position = 'none')

##############################################################################



#3
###################################################################################
################### data group by histological subtype ############################

##(1)
### turn into function()
# only need tumour_id
generate_histological_group <- function(tumortype){
  histological_tumour_id = tumor_df[ Histology_per_tumour_id_muttable==tumortype][,c("tumour_id")]
  histological_tumour_id <-as.data.table( histological_tumour_id)[,c("tumour_id")]
  return(histological_tumour_id)
}

##(2)
###  LOOP to run all subtypes into a mutation list
#for loop to make all subtype mutation
subtype <- c("Invasive adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma",
             "Pleomorphic carcinoma","LCNEC","Large cell carcinoma","combined LUAD and LCNEC",
             "Collision LUAD and LUSC","Carcinosarcoma")
#Note the opening case
                                                                              
histological_mutation_list <- list()
for (x in subtype){
  table <- generate_histological_group(x )
  histological_mutation_list[[x]] <- table
}

head(histological_mutation_list)

#############################################################################




#4
###################################################################################
################### make mutation tables grouped by tumor type ##################
#(1)
######### LOOP to merge all subtype tables with corresponding mutations
for (i in 1:length(histological_mutation_list)) {
  table <- histological_mutation_list[[i]]
  merged_table <- merge(table, mutation, by = "tumour_id")
  histological_mutation_list[[i]] <- merged_table
}

head(histological_mutation_list)

##################################################################################




#5
###################################################################################
######################### run dndscv on mutation tables ###########################
#(1)
##### LOOP to run dndscv and output all dndsout_subtype 
short_name <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")
subtype <- c("Invasive adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma",
             "Pleomorphic carcinoma","LCNEC","Large cell carcinoma","combined LUAD and LCNEC",
             "Collision LUAD and LUSC","Carcinosarcoma")
library(dndscv)

dndsout_list <- list()
for (i in 1:length(subtype)){
    sublist <- histological_mutation_list [[i]]
    list <- dndscv(sublist, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)
    dndsout_list[[i]] <- assign(paste0("dndsout_", short_name[i]), list, envir = .GlobalEnv)
}  

## Warning messages:(for Collision LUAD and LUSC) In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  : iteration limit reached


#(2)
##### LOOP to output significant genes
subtype <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")

for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(head(sel_cv_subtype), digits = 3) 
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}  


#(3) 
###### LOOP to get Global dN/dS estimates[globaldnds] for every subtype

for (i in subtype) {
  print(i)
  print("Global dN/dS")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  print (dndsout_subtype$globaldnds)
}  


#(4)
##### other outputs
# 【annotmuts】an annotated table of coding mutations
# 【mle_submodel】MLEs of mutation rate parameters
# 【genemuts】 a table with the observed expected number of mutations per gene 

#(5)
##### LOOP to get theta/θ value to test suitability for dNdScv

for (i in subtype) {
  print(i)
  print("theta value")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(dndsout_subtype$nbreg$theta)
}

#(6)
#####  local neutrality test

for (i in subtype) {
  print(i)
  print("local neutrality test")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  signif_genes_localmodel_subtype = as.vector(dndsout_subtype$sel_loc$gene_name[dndsout_subtype$sel_loc$qall_loc<0.1])
  print(signif_genes_localmodel_subtype)
}


#(7)
##### AIC model to measure the fit of a statistical model 
# the smaller the AIC, the better the model

for (i in subtype) {
  print(i)
  print("AIC")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(AIC(dndsout_subtype$poissmodel))
}

#############################################################







