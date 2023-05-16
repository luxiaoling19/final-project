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



#3
###################################################################################
################### data group by histological subtype ############################
#Squamous_cell_carcinoma_tumours= tumor_df[ Histology_per_tumour_id_muttable=="Squamous cell carcinoma"]
#Squamous_cell_carcinoma_tumours <- as.data.table(Squamous_cell_carcinoma_tumours)[,c("tumour_id")]

##(1)
### turn into function()
generate_histological_group <- function(tumortype){
  histological_tumour_id = tumor_df[ Histology_per_tumour_id_muttable==tumortype][,c("tumour_id")]
  histological_tumour_id <-as.data.table( histological_tumour_id)[,c("tumour_id")]
  return(histological_tumour_id)
}
# 调用function：input histological tumor type, out put _tumor_id
#LCNEC_tumour_id <- generate_histological_group("LCNEC")
# only need tumour_id

##(2)
##########  LOOP to run all subtypes into a mutation list
#for loop to make all subtype mutation
subtype <- c("Invasive adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma",
             "Pleomorphic carcinoma","LCNEC","Large cell carcinoma","combined LUAD and LCNEC",
             "Collision LUAD and LUSC","Carcinosarcoma")
#Note the opening case注意开头大小写
                                                                              
histological_mutation_list <- list()
for (x in subtype){
  table <- generate_histological_group(x )
  histological_mutation_list[[x]] <- table
}

head(histological_mutation_list)

#LCNEC_mutation <- histological_mutation_list$LCNEC

#############################################################################




#4
###################################################################################
################### make mutation tables grouped by tumor type ##################
# Adenosquamous_carcinoma_mutation <- merge(Adenosquamous_carcinoma_tumours,mutation,by="tumour_id")

#turn into function()
# merge_histological_type_with_mutation <- function(histological_table){merge(histological_table, mutation,by="tumour_id")}
#implement function()
# Squamous_cell_carcinoma_mutation <- merge_histological_type_with_mutation(Squamous_cell_carcinoma_tumours)


######### LOOP to merge all subtype tables with corresponding mutations
for (i in 1:length(histological_mutation_list)) {
  table <- histological_mutation_list[[i]]
  merged_table <- merge(table, mutation, by = "tumour_id")
  histological_mutation_list[[i]] <- merged_table
}

head(histological_mutation_list)

##################################################################################



############# remove data and values do not need ###############################
objects_to_remove <- c("patient_by_subtype","merged_table","x","i","table","mutation","tumor_df","subtype")
rm(list = objects_to_remove, envir = .GlobalEnv)
# 检查数据对象和变量是否成功删除
ls()



#5
###################################################################################
######################### run dndscv on mutation tables ###########################

########## run dndscv
# Adenosquamous_carcinoma_mutation$chr = gsub("chr","",as.vector(Adenosquamous_carcinoma_mutation$chr))

##!!! for invasive adenocarcinoma--infinity
#Note: 1 samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). 247 samples left after filtering.
#Note: 458 mutations removed for exceeding the limit of mutations per gene per sample (see the max_muts_per_gene_per_sample argument in dndscv)
## !!sample too big,should set "max_muts_per_gene_per_sample" and "max_coding_muts_per_sample" to infinity/Inf)

#dndsout_IAC = dndscv(histological_mutation_list$`Invasive adenocarcinoma`, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)


############### LOOP to run dndscv and output all dndsout_subtype #################
short_name <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")
subtype <- c("Invasive adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma",
             "Pleomorphic carcinoma","LCNEC","Large cell carcinoma","combined LUAD and LCNEC",
             "Collision LUAD and LUSC","Carcinosarcoma")


dndsout_list <- list()
for (i in 1:length(subtype)){
    sublist <- histological_mutation_list [[i]]
    list <- dndscv(sublist, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)
    dndsout_list[[i]] <- assign(paste0("dndsout_", short_name[i]), list, envir = .GlobalEnv)

}  




##!!! for Collision LUAD and LUSC
# ????? Warning messages:In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  : iteration limit reached



############### LOOP to get significant genes for every subtype #################

subtype <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")

#(1) 
############# output significant genes 输出显著性基因 ##############
#LOOP to get significant genes for every subtype
for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(head(sel_cv_subtype), digits = 3) 
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}  


for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv","wmis_cv","wnon_cv","wspl_cv","wind_cv")], digits = 3)
}  


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

#(5)
############ 局部中性测试 local neutrality test #############
#LOOP to get local neutrality test for every subtype
for (i in subtype) {
  print(i)
  print("local neutrality test")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  signif_genes_localmodel_subtype = as.vector(dndsout_subtype$sel_loc$gene_name[dndsout_subtype$sel_loc$qall_loc<0.1])
  print(signif_genes_localmodel_subtype)
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


###############################
#high MLE of the dN/dS ratio + significant gene = higher chance of driver mutation






