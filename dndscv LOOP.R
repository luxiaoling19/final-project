#(1) 
######## output significant genes 输出显著性基因 #####
sel_cv_ASC = dndsout_ASC $ sel_cv
print(head(sel_cv_ASC), digits = 3) 

signif_genes_ASC = sel_cv_ASC[sel_cv_ASC $ qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
head(signif_genes_ASC)

# 使用geneci计算置信区间，默认Confidence level = 0.95
dndsout_ASC2 <- dndscv(histological_mutation_list$`Adenosquamous carcinoma`,outmats=T)
ci_ASC2 = geneci(dndsout_ASC2)
#geneci(dndsout, gene_list = NULL, level = 0.95)

#(2) 
########## 输出全局dN/dS Global dN/dS estimates [写作globaldnds] ############
print (dndsout_ASC$globaldnds)
#【wmis】missense,【wnon】nonsense,【wspl】essential splice site substitutions【wall】all non-synonymous substitutions
#【wtru】truncating substitutions together(include nonsense and essential splice site mutations)
# Normally, all should be ~1. If dN/dS<<1， reflected a problem of SNP contamination or an inadequate substitution model

#(3)
######## other outputs ########
head(dndsout_ASC$annotmuts)
# 【annotmuts】an annotated table of coding mutations编码突变的注释表 
# 【mle_submodel】MLEs of mutation rate parameters突变率参数的最大似然估计 
# 【genemuts】 a table with the observed expected number of mutations per gene 每个基因观察到的和预期的突变数量表

#(4)
######### theta/θ value to test suitability for dNdScv #######
# theta_ASC <- dndsout_ASC$nbreg$theta
print(dndsout_ASC$nbreg$theta)
##Result：very low θ estimate(过分散参数the overdispersion parameter)，especially θ < 1，dNdScv not applicable

#(5)
############ 局部中性测试 local neutrality test #############
signif_genes_localmodel_ASC = as.vector(dndsout_ASC$sel_loc$gene_name[dndsout_ASC$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel_ASC)

#(6)
############### AIC model to measure the fit of a statistical model #############
# the smaller the AIC, the better the model
AIC(dndsout_ASC$poissmodel)



############### LOOP to get significant genes for every subtype #################
subtype <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")

for (i in subtype) {
  print(i)
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}


subtype <- c("IAC", "SCC", "ASC", "PPC", "LCNEC", "LCC", "combined", "Collision", "CS")
subtype <-c("ASC")
#LOOP to get significant genes for every subtype
for (i in subtype) {
  print(i)
  print("significant genes")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  sel_cv_subtype <- dndsout_subtype$sel_cv
  print(sel_cv_subtype[sel_cv_subtype$qglobal_cv < 0.1, c("gene_name", "qglobal_cv")], digits = 3)
}  
#LOOP to get Global dN/dS for every subtype
for (i in subtype) {
  print(i)
  print("Global dN/dS")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  print (dndsout_subtype$globaldnds)
}  
#LOOP to get theta value for every subtype
for (i in subtype) {
  print(i)
  print("theta value")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(dndsout_subtype$nbreg$theta)
}
#LOOP to get local neutrality test for every subtype
for (i in subtype) {
  print(i)
  print("local neutrality test")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)
  signif_genes_localmodel_subtype = as.vector(dndsout_subtype$sel_loc$gene_name[dndsout_subtype$sel_loc$qall_loc<0.1])
  print(signif_genes_localmodel_ASC)
}
#LOOP to get AIC for every subtype
for (i in subtype) {
  print(i)
  print("AIC")
  table <- paste0("dndsout_", i) 
  dndsout_subtype <- get(table)  
  print(AIC(dndsout_subtype$poissmodel))
}

