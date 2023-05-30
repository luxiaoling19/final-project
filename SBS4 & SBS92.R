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
################  SBS4>0.1 and SBS92 >0.1 ##########################

#(1)
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
rm(signature_weights)


signature_weights_perTumour <- signature_weights_perTumour[signature_weights_perTumour[, "SBS4"] > 0.1 & signature_weights_perTumour[, "SBS92"] > 0.1, 
                                                           "tumour_id", drop = FALSE] 
#drop = FALSE:Output result in the form of a data frame (not a vector)
head(signature_weights_perTumour)

################################################################################
SBS4_SBS92_both_0.1 <- merge(signature_weights_perTumour,mutation, by = "tumour_id")
head(SBS4_SBS92_both_0.1)
################################################################################
dndsout_SBS4_SBS92_both_0.1 <- dndscv(SBS4_SBS92_both_0.1, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)

sel_cv = dndsout_SBS4_SBS92_both_0.1 $ sel_cv
signif_genes = sel_cv[sel_cv $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
print(signif_genes,digit=3)

print (dndsout_SBS4_SBS92_both_0.1$globaldnds)

print(dndsout_SBS4_SBS92_both_0.1$nbreg$theta)

AIC(dndsout_SBS4_SBS92_both_0.1$poissmodel)
################################################################################






#2
################################################################################
################  SBS4>0.1 OR SBS92 >0.1 ##########################
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
rm(signature_weights)

SBS4_plus_SBS92 <- signature_weights_perTumour[signature_weights_perTumour[, "SBS4"] > 0.1 | signature_weights_perTumour[, "SBS92"] > 0.1,"tumour_id",drop = FALSE]
head(SBS4_plus_SBS92)
###################################################################
SBS4_plus_SBS92 <- merge(SBS4_plus_SBS92,mutation, by = "tumour_id")
head(SBS4_plus_SBS92)
###################################################################
dndsout_SBS4_plus_SBS92 <- dndscv(SBS4_plus_SBS92, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)

sel_cv = dndsout_SBS4_plus_SBS92 $ sel_cv
signif_genes = sel_cv[sel_cv $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
print(signif_genes,digit=3)

print (dndsout_SBS4_plus_SBS92$globaldnds)

print(dndsout_SBS4_plus_SBS92$nbreg$theta)

AIC(dndsout_SBS4_plus_SBS92$poissmodel)