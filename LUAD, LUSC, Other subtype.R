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

#subtype <- c("Invasive adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma", "Pleomorphic carcinoma","LCNEC","Large cell carcinoma","combined LUAD and LCNEC", "Collision LUAD and LUSC","Carcinosarcoma")

##############################################################################



#3 group patients by LUAD, LUSC, other types
###################################################################################
################### data group by histological subtype ############################


##(1)
### turn into function()
generate_histological_group <- function(tumortype){
  histological_tumour_id = tumor_df[Histology_per_tumour_id_muttable %in% tumortype][,c("tumour_id")]
  histological_tumour_id <-as.data.table( histological_tumour_id)[,c("tumour_id")]
  return(histological_tumour_id)
}

##(2)
##########  LOOP to run all subtypes into a mutation list

subtype1 <- "Invasive adenocarcinoma"
subtype2 <- "Squamous cell carcinoma"
other_subtypes <- c("Adenosquamous carcinoma","Pleomorphic carcinoma",
             "LCNEC","Large cell carcinoma","combined LUAD and LCNEC",
             "Collision LUAD and LUSC","Carcinosarcoma")

histological_mutation_list <- list()
  
table <- generate_histological_group(subtype1)
histological_mutation_list[["LUAD"]] <- table
table <- generate_histological_group(subtype2)
histological_mutation_list[["LUSC"]] <- table 
table <- generate_histological_group(other_subtypes)
histological_mutation_list[["Other subtypes"]] <- table

head(histological_mutation_list)

#############################################################################




#4
###################################################################################
################### make mutation tables grouped by tumor type ##################


######### LOOP to merge all subtype tables with corresponding mutations
for (i in 1:length(histological_mutation_list)) {
  table <- histological_mutation_list[[i]]
  merged_table <- merge(table, mutation, by = "tumour_id")
  histological_mutation_list[[i]] <- merged_table
}

head(histological_mutation_list)

##################################################################################



############# remove data and values do not need ###############################
# objects_to_remove <- c("patient_by_subtype","merged_table","x","i","table","mutation","tumor_df","subtype")
# rm(list = objects_to_remove, envir = .GlobalEnv)
# ls()



#5
###################################################################################
######################### run dndscv on mutation tables ###########################

subtype <- c("LUAD","LUSC","Other subtypes")
dndsout_list <- list()
for (i in 1:length(subtype)){
  sublist <- histological_mutation_list [[i]]
  list <- dndscv(sublist, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)
  dndsout_list[[i]] <- assign(paste0("dndsout_", subtype[i]), list, envir = .GlobalEnv)
}  


####################### LOOP to get output #########################

#(1) 
############# output significant genes 输出显著性基因 ##############
#LOOP to get significant genes for every subtype

histological_driver_q0.1_list <- list()
for (i in subtype ) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  histological_driver_q0.1_list[[i]] <- sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")]
}  

#save driver genes to excel document
library(writexl)
folder_path <- "D:/Program Files/R/final project"
# 确保文件夹存在，如果不存在则创建
#dir.create(folder_path, showWarnings = FALSE)
excel_file <- file.path(folder_path, "histological_driver_q0.1_list.xlsx")
write_xlsx(histological_driver_q0.1_list, excel_file)


#(2) 
########## 输出全局dN/dS Global dN/dS estimates [写作globaldnds] ############
#LOOP to get Global dN/dS for every subtype
for (i in subtype) {
  print(i)
  print("Global dN/dS")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  print (dndsout_subtype$globaldnds)
}  


#(3)
######## other outputs ########
head(dndsout_ASC$annotmuts)
# 【annotmuts】an annotated table of coding mutations编码突变的注释表 
# 【mle_submodel】MLEs of mutation rate parameters突变率参数的最大似然估计 
# 【genemuts】 a table with the observed expected number of mutations per gene 每个基因观察到的和预期的突变数量表
annotmuts_list <-  list()
for (i in subtype) {
  print(i)
  print("annotated table of coding mutations")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  annotmuts_list[[table]] <- dndsout_subtype$annotmuts
}
observed_expected_number_of_mutations <- list()
for (i in subtype) {
  print(i)
  print("a table with the observed expected number of mutations per gene")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  observed_expected_number_of_mutations[[table]] <- dndsout_subtype$genemuts
}

#(4)
######### theta/θ value to test suitability for dNdScv #######
#LOOP to get theta value for every subtype
for (i in subtype) {
  print(i)
  print("theta value")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(dndsout_subtype$nbreg$theta)
}


#(6)
############### AIC model to measure the fit of a statistical model #############
# the smaller the AIC, the better the model
#LOOP to get AIC for every subtype
for (i in subtype) {
  print(i)
  print("AIC")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(AIC(dndsout_subtype$poissmodel))
}

########################################################################
