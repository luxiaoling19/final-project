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
################  SBS2>0.1 and SBS13 >0.1 ##########################

#(1)
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
rm(signature_weights)


signature_weights_perTumour <- signature_weights_perTumour[signature_weights_perTumour[, "SBS2"] > 0.1 & signature_weights_perTumour[, "SBS13"] > 0.1, 
                                                           "tumour_id", drop = FALSE] 
#drop = FALSE:Output result in the form of a data frame (not a vector)
head(signature_weights_perTumour)

################################################################################
SBS2_SBS13_both_0.1 <- merge(signature_weights_perTumour,mutation, by = "tumour_id")
head(SBS2_SBS13_both_0.1)
################################################################################
dndsout_SBS2_SBS13_both_0.1 <- dndscv(SBS2_SBS13_both_0.1, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)

sel_cv = dndsout_SBS2_SBS13_both_0.1 $ sel_cv
signif_genes = sel_cv[sel_cv $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
print(signif_genes,digit=3)

print (dndsout_SBS2_SBS13_both_0.1$globaldnds)

print(dndsout_SBS2_SBS13_both_0.1$nbreg$theta)

AIC(dndsout_SBS2_SBS13_both_0.1$poissmodel)
################################################################################






#2
################################################################################
################  SBS2>0.1 OR SBS13 >0.1 ##########################
signature_weights = readRDS('figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds')
signature_weights_perTumour <- signature_weights[["signature_weights_perTumour"]]
rm(signature_weights)

SBS2_plus_SBS13 <- signature_weights_perTumour[signature_weights_perTumour[, "SBS2"] > 0.1 | signature_weights_perTumour[, "SBS13"] > 0.1, 
                                                           "tumour_id", drop = FALSE]
head(SBS2_plus_SBS13)
###################################################################
SBS2_plus_SBS13 <- merge(SBS2_plus_SBS13,mutation, by = "tumour_id")
head(SBS2_plus_SBS13)
###################################################################
dndsout_SBS2_plus_SBS13 <- dndscv(SBS2_plus_SBS13, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)

sel_cv = dndsout_SBS2_plus_SBS13 $ sel_cv
signif_genes = sel_cv[sel_cv $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
print(signif_genes,digit=3)

print (dndsout_SBS2_plus_SBS13$globaldnds)

print(dndsout_SBS2_plus_SBS13$nbreg$theta)

AIC(dndsout_SBS2_plus_SBS13$poissmodel)



###################################################################
major_SBS2_OR_SBS13 <- max_SBS[grepl("SBS2|SBS13",max_SBS$max_sig_weights),"tumour_id"]

major_SBS2_or_SBS13 <- merge(major_SBS2_OR_SBS13,mutation, by = "tumour_id")
head(major_SBS2_or_SBS13)
dndsout_major_SBS2_or_SBS13 <- dndscv(major_SBS2_or_SBS13, max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf)

sel_cv = dndsout_major_SBS2_or_SBS13 $ sel_cv
signif_genes = sel_cv[sel_cv $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
print(signif_genes,digit=3)

print (dndsout_major_SBS2_or_SBS13$globaldnds)

print(dndsout_major_SBS2_or_SBS13$nbreg$theta)

AIC(dndsout_major_SBS2_or_SBS13$poissmodel)
