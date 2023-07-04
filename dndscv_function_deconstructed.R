#!/usr/bin/env Rscript
# FILE: dndscv_function_deconstructed.R ---------------------------------------
#
# DESCRIPTION: this script contains dndscv function from dNdScv package 
#              annotated by Maria Litovchenko
# USAGE: Read the script and run it line by line, try to undertand what it does
#        and how it does it
#
# OPTIONS: 
#
# EXAMPLE:
#
# REQUIREMENTS: 
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  1
# CREATED:  04.07.2023
# REVISION: 04.07.2023

library(dndscv)
# Input arguments -------------------------------------------------------------
# In this section input test arguments for dndscv function are listed. They 
# are here so you could run the script and see what it does.

# to load standard mutations table which comes together with the package
data("dataset_simbreast", package="dndscv")

gene_list = NULL
refdb = "hg19"
sm = "192r_3w"
kc = "cgc81"
cv = "hg19"
max_muts_per_gene_per_sample = 3
max_coding_muts_per_sample = 3000
use_indel_sites = T
min_indels = 5
maxcovs = 20
constrain_wnon_wspl = T
outp = 3
numcode = 1 
outmats = F
mingenecovs = 500

# Functions -------------------------------------------------------------------
# All this functions originally were written inside dnsdcv function. Yes, you 
# can write a function inside the function. But IT IS REALLY BAD PRACTICE!!!
# DO NOT DO THAT! Therefore, I moved them here, separetly as all little good
# functions should be

chr2cds <-  function(pos, cds_int, strand) {
  if (strand == 1) {
    return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
  }
  else if (strand == -1) {
    return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% 
                   pos))
  }
}

fit_substmodel = function(N, L, substmodel) {
  l = c(L)
  n = c(N)
  r = c(substmodel)
  n = n[l != 0]
  r = r[l != 0]
  l = l[l != 0]
  params = unique(base::strsplit(x = paste(r, collapse = "*"), 
                                 split = "\\*")[[1]])
  indmat = as.data.frame(array(0, dim = c(length(r), length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) {
    indmat[j, base::strsplit(r[j], split = "\\*")[[1]]] = 1
  }
  model = glm(formula = n ~ offset(log(l)) + . - 1, data = indmat, 
              family = poisson(link = log))
  mle = exp(coefficients(model))
  ci = exp(confint.default(model))
  par = data.frame(name = gsub("`", "", rownames(ci)), mle = mle[rownames(ci)],
                   cilow = ci[, 1], cihigh = ci[, 2])
  return(list(par = par, model = model))
}

selfun_loc = function(j) {
  x = RefCDS[[j]]
  y = as.numeric(genemuts[j, -1])
  mrfold = sum(y[1:4])/sum(y[5:8])
  ll0 = sum(dpois(x = x$N, lambda = x$L * mutrates * 
                    mrfold * t(array(c(1, 1, 1, 1), dim = c(4, numrates))), 
                  log = T))
  mrfold = max(1e-10, sum(y[c(1, 2)])/sum(y[c(5, 6)]))
  wfree = y[3:4]/y[7:8]/mrfold
  wfree[y[3:4] == 0] = 0
  llmis = sum(dpois(x = x$N, 
                    lambda = x$L * mutrates *  mrfold * 
                      t(array(c(1, 1, wfree), dim = c(4, numrates))), 
                    log = T))
  mrfold = max(1e-10, y[1]/y[5])
  w = y[2:4]/y[6:8]/mrfold
  w[y[2:4] == 0] = 0
  llall = sum(dpois(x = x$N, 
                    lambda = x$L * mutrates * mrfold * t(array(c(1, w), 
                                                               dim = c(4, numrates))), 
                    log = T))
  w[w > 10000] = 10000
  p = 1 - pchisq(2 * (llall - c(llmis, ll0)), df = c(1, 3))
  return(c(w, p))
}

mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
  tml = (n_neutral + shape - 1)/(exp_rel_neutral + (1/scale))
  if (shape <= 1) {
    tml = max(shape * scale, tml)
  }
  return(tml)
}

selfun_cv = function(j) {
  x = RefCDS[[j]]
  y = as.numeric(genemuts[j, -1])
  exp_rel = y[5:8]/y[5]
  shape = theta
  scale = y[9]/theta
  indneut = 1:4
  opt_t = mle_tcv(n_neutral = sum(y[indneut]), shape = shape, scale = scale,
                  exp_rel_neutral = sum(exp_rel[indneut]))
  mrfold = max(1e-10, opt_t/y[5])
  ll0 = sum(dpois(x = x$N, lambda = x$L * mutrates * 
                    mrfold * t(array(c(1, 1, 1, 1), dim = c(4, numrates))), 
                  log = T)) + dgamma(opt_t, shape = shape, scale = scale, 
                                     log = T)
  indneut = 1:2
  opt_t = mle_tcv(n_neutral = sum(y[indneut]), shape = shape, scale = scale,
                  exp_rel_neutral = sum(exp_rel[indneut]))
  mrfold = max(1e-10, opt_t/sum(y[5]))
  wfree = y[3:4]/y[7:8]/mrfold
  wfree[y[3:4] == 0] = 0
  llmis = sum(dpois(x = x$N, lambda = x$L * mutrates * mrfold * t(array(c(1, 1, wfree), dim = c(4, numrates))), 
                    log = T)) + dgamma(opt_t, shape = shape, scale = scale, 
                                       log = T)
  indneut = c(1, 3, 4)
  opt_t = mle_tcv(n_neutral = sum(y[indneut]), shape = shape, scale = scale,
                  exp_rel_neutral = sum(exp_rel[indneut]))
  mrfold = max(1e-10, opt_t/sum(y[5]))
  wfree = y[2]/y[6]/mrfold
  wfree[y[2] == 0] = 0
  lltrunc = sum(dpois(x = x$N, lambda = x$L * mutrates * 
                        mrfold * t(array(c(1, wfree, 1, 1), dim = c(4, 
                                                                    numrates))), log = T)) + dgamma(opt_t, shape = shape, 
                                                                                                    scale = scale, log = T)
  indneut = 1
  opt_t = mle_tcv(n_neutral = sum(y[indneut]), shape = shape, scale = scale,
                  exp_rel_neutral = sum(exp_rel[indneut]))
  mrfold = max(1e-10, opt_t/sum(y[5]))
  wfree = y[2:4]/y[6:8]/mrfold
  wfree[y[2:4] == 0] = 0
  llall_unc = sum(dpois(x = x$N, lambda = x$L * mutrates * 
                          mrfold * t(array(c(1, wfree), dim = c(4, numrates))), 
                        log = T)) + dgamma(opt_t, shape = shape, scale = scale, 
                                           log = T)
  if (constrain_wnon_wspl == 0) {
    p = 1 - pchisq(2 * (llall_unc - c(llmis, lltrunc, ll0)), df = c(1, 2, 3))
    return(c(wfree, p))
  }
  else {
    wmisfree = y[2]/y[6]/mrfold
    wmisfree[y[2] == 0] = 0
    wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold
    wtruncfree[sum(y[3:4]) == 0] = 0
    llall = sum(dpois(x = x$N, lambda = x$L * mutrates * 
                        mrfold * t(array(c(1, wmisfree, wtruncfree, 
                                           wtruncfree), dim = c(4, numrates))), log = T)) + 
      dgamma(opt_t, shape = shape, scale = scale, 
             log = T)
    p = 1 - pchisq(2 * c(llall_unc - llmis, llall - 
                           c(lltrunc, ll0)), df = c(1, 1, 2))
    return(c(wmisfree, wtruncfree, wtruncfree, p))
  }
}

# Read mutations --------------------------------------------------------------
message("[1] Loading the environment...")
mutations = mutations[, 1:5]
mutations[, c(1, 2, 3, 4, 5)] = lapply(mutations[, c(1, 2, 3, 4, 5)],
                                       as.character)
mutations[[3]] = as.numeric(mutations[[3]])
mutations = mutations[mutations[, 4] != mutations[, 5], ]
colnames(mutations) = c("sampleID", "chr", "pos", "ref", "mut")
indna = which(is.na(mutations), arr.ind = T)
if (nrow(indna) > 0) {
  mutations = mutations[-unique(indna[, 1]), ]
  warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.", 
                  length(unique(indna[, 1]))))
}

# Read positions of genes -----------------------------------------------------
refdb_class = class(refdb)
if ("character" %in% refdb_class) {
  if (refdb == "hg19") {
    data("refcds_hg19", package = "dndscv")
    if (any(gene_list == "CDKN2A")) {
      gene_list = unique(c(setdiff(gene_list, "CDKN2A"), 
                           "CDKN2A.p14arf", "CDKN2A.p16INK4a"))
    }
  }
  else {
    load(refdb)
  }
} else if ("array" %in% refdb_class) {
  RefCDS = refdb
} else {
  stop("Expected refdb to be \"hg19\", a file path, or a RefCDS-formatted array object.")
}

# SUPER IMPORTANT NOTE --------------------------------------------------------
# RefCDS is essentially a list. Each element of a list contains informations 
# about one gene, including its name, gene and protein ids, strand, sequence of
# CDS, etc. Also, it contains (for every gene) a matrix called L. This is a 
# matrix with 192 rows, each corresponding to a possible 3-nucleotide change
# and 4 columns. Each column corresponds to a one of possible impacts of 
# mutation on CDS: first column - synonymous, second column - nonsense, 
# third column - missense, forth - anything else. The value in each cell of 
# matrix corresponds to the number of possible mutations with such impact in
# the certain trinucleotide context across the whole gene.

# Read gene list --------------------------------------------------------------
if (is.null(gene_list)) {
  gene_list = sapply(RefCDS, function(x) x$gene_name)
} else {
  allg = sapply(RefCDS, function(x) x$gene_name)
  nonex = gene_list[!(gene_list %in% allg)]
  if (length(nonex) > 0) {
    stop(sprintf("The following input gene names are not in the RefCDS database: %s", 
                 paste(nonex, collapse = ", ")))
  }
  RefCDS = RefCDS[allg %in% gene_list]
  gr_genes = gr_genes[gr_genes$names %in% gene_list]
}

# Read covariates -------------------------------------------------------------
if (is.character(cv)) {
  data(list = sprintf("covariates_%s", cv), package = "dndscv")
} else {
  covs = cv
}

# Read cancer genes -----------------------------------------------------------
if (kc[1] %in% c("cgc81")) {
  data(list = sprintf("cancergenes_%s", kc), package = "dndscv")
} else {
  known_cancergenes = kc
}

# Read substitution model -----------------------------------------------------
if (length(sm) == 1) {
  data(list = sprintf("submod_%s", sm), package = "dndscv")
} else {
  substmodel = sm
}

# Convert CDS to array of letters ---------------------------------------------
for (j in 1:length(RefCDS)) {
  # converts CDS from a string to array of individual letters
  RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), 
                                       split = "")[[1]]
  # converts CDS shifted up by 1 letter from a string to array of individual
  # letters
  RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), 
                                          split = "")[[1]]
  # converts CDS shifted down by 1 letter from a string to array of individual
  # letters
  RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), 
                                            split = "")[[1]]
  # do the same for splice sites
  if (!is.null(RefCDS[[j]]$seq_splice)) {
    RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), 
                                            split = "")[[1]]
    RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), 
                                               split = "")[[1]]
    RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), 
                                                 split = "")[[1]]
  }
}

# Convert mutations for GRanges object ----------------------------------------
message("[2] Annotating the mutations...")
# assign an index for each gene
ind = setNames(1:length(RefCDS), sapply(RefCDS, function(x) x$gene_name))
gr_genes_ind = ind[gr_genes$names]
if (any(diff(mutations$pos) == 1)) {
  warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
}
if (nrow(unique(mutations[, 2:5])) < nrow(mutations)) {
  warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
}
mutations$end = mutations$start = mutations$pos
# length of mutation
l = nchar(mutations$ref) - 1
mutations$end = mutations$end + l
# this vector will hold true, if mutation is indel, and false otherwise
ind = substr(mutations$ref, 1, 1) == substr(mutations$mut, 1, 1) & 
  nchar(mutations$ref) > nchar(mutations$mut)
mutations$start = mutations$start + ind
gr_muts = GenomicRanges::GRanges(mutations$chr, 
                                 IRanges::IRanges(mutations$start,
                                                  mutations$end))

# Annotate mutations to genes -------------------------------------------------
ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type = "any",
                                               select = "all"))
mutations = mutations[ol[, 1], ]
mutations$geneind = gr_genes_ind[ol[, 2]]
mutations$gene = sapply(RefCDS, function(x) x$gene_name)[mutations$geneind]
mutations = unique(mutations)

# count number of mutations per sample & remove samples which are too 
# high mutated
nsampl = sort(table(mutations$sampleID))
exclsamples = NULL
if (any(nsampl > max_coding_muts_per_sample)) {
  message(sprintf("    Note: %0.0f samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). %0.0f samples left after filtering.", 
                  sum(nsampl > max_coding_muts_per_sample), 
                  sum(nsampl <= max_coding_muts_per_sample)))
  exclsamples = names(nsampl[nsampl > max_coding_muts_per_sample])
  mutations = mutations[!(mutations$sampleID %in% 
                            names(nsampl[nsampl > max_coding_muts_per_sample])), ]
}
mutrank = ave(mutations$pos, paste(mutations$sampleID, mutations$gene), 
              FUN = function(x) rank(x))
exclmuts = NULL
if (any(mutrank > max_muts_per_gene_per_sample)) {
  message(sprintf("    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample (see the max_muts_per_gene_per_sample argument in dndscv)", 
                  sum(mutrank > max_muts_per_gene_per_sample)))
  exclmuts = mutations[mutrank > max_muts_per_gene_per_sample, ]
  mutations = mutations[mutrank <= max_muts_per_gene_per_sample, ]
}
mutations$strand = sapply(RefCDS, function(x) x$strand)[mutations$geneind]

# Split SNVs from indels ------------------------------------------------------
nt = c("A", "C", "G", "T")
# complimentary nucleotides
compnt = setNames(rev(nt), nt)

snv = (mutations$ref %in% nt & mutations$mut %in% nt)
if (!any(snv)) {
  stop("Zero coding substitutions found in this dataset. Unable to run dndscv. Common causes for this error are inputting only indels or using chromosome names different to those in the reference database (e.g. chr1 vs 1)")
}
indels = mutations[!snv, ]
mutations = mutations[snv, ]
mutations$ref_cod = mutations$ref
mutations$mut_cod = mutations$mut
# for the mutations on the opposite strand, muts as reference as mutated 
# nucleotides there compliments
isminus = (mutations$strand == -1)
mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]

# Create array which will hold all possible 192 3-nucleotide changes ----------
nt = c("A", "C", "G", "T")
trinucs = paste(rep(nt, each = 16, times = 1), rep(nt, each = 4,  times = 4), 
                rep(nt, each = 1, times = 16), sep = "")
trinucinds = setNames(1:64, trinucs)
trinucsubs = NULL
for (j in 1:length(trinucs)) {
  trinucsubs = c(trinucsubs, 
                 paste(trinucs[j], 
                       paste(substr(trinucs[j], 1, 1), 
                             setdiff(nt, substr(trinucs[j], 2, 2)),
                             substr(trinucs[j], 3, 3), sep = ""), sep = ">"))
}
trinucsubsind = setNames(1:192, trinucsubs)

# Annotate SNVs with their change of coding region ----------------------------
for (j in 1:length(RefCDS)) {
  RefCDS[[j]]$N = array(0, dim = c(192, 4))
}
ref3_cod = mut3_cod = wrong_ref = aachange = array(NA, nrow(mutations))
ntchange = impact = codonsub = array(NA, nrow(mutations))

for (j in 1:nrow(mutations)) {
  geneind = mutations$geneind[j]
  pos = mutations$pos[j]
  if (any(pos == RefCDS[[geneind]]$intervals_splice)) {
    impact[j] = "Essential_Splice"
    impind = 4
    pos_ind = (pos == RefCDS[[geneind]]$intervals_splice)
    cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
    ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], 
                          RefCDS[[geneind]]$seq_splice[pos_ind], 
                          RefCDS[[geneind]]$seq_splice1down[pos_ind])
    mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], 
                          mutations$mut_cod[j],
                          RefCDS[[geneind]]$seq_splice1down[pos_ind])
    aachange[j] = ntchange[j] = codonsub[j] = "."
  }
  else {
    pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, 
                      RefCDS[[geneind]]$strand)
    cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
    ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], 
                          RefCDS[[geneind]]$seq_cds[pos_ind], 
                          RefCDS[[geneind]]$seq_cds1down[pos_ind])
    mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], 
                          mutations$mut_cod[j], 
                          RefCDS[[geneind]]$seq_cds1down[pos_ind])
    codon_pos = c(ceiling(pos_ind/3) * 3 - 2, ceiling(pos_ind/3) * 3 - 1, 
                  ceiling(pos_ind/3) * 3)
    old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
    pos_in_codon = pos_ind - (ceiling(pos_ind/3) - 1) * 3
    new_codon = old_codon
    new_codon[pos_in_codon] = mutations$mut_cod[j]
    old_aa = seqinr::translate(old_codon, numcode = numcode)
    new_aa = seqinr::translate(new_codon, numcode = numcode)
    aachange[j] = sprintf("%s%0.0f%s", old_aa, ceiling(pos_ind/3), new_aa)
    ntchange[j] = sprintf("%s%0.0f%s", mutations$ref_cod[j], pos_ind, 
                          mutations$mut_cod[j])
    codonsub[j] = sprintf("%s>%s", paste(old_codon, collapse = ""), 
                          paste(new_codon, collapse = ""))
    if (new_aa == old_aa) {
      impact[j] = "Synonymous"
      impind = 1
    }
    else if (new_aa == "*") {
      impact[j] = "Nonsense"
      impind = 3
    }
    else if (old_aa != "*") {
      impact[j] = "Missense"
      impind = 2
    }
    else if (old_aa == "*") {
      impact[j] = "Stop_loss"
      impind = NA
    }
  }
  if (mutations$ref_cod[j] != as.character(cdsnt)) {
    wrong_ref[j] = 1
  }
  else if (!is.na(impind)) {
    trisub = trinucsubsind[paste(ref3_cod[j], mut3_cod[j], sep = ">")]
    RefCDS[[geneind]]$N[trisub, impind] = RefCDS[[geneind]]$N[trisub, 
                                                              impind] + 1
  }
  if (round(j/10000) == (j/10000)) {
    message(sprintf("    %0.3g%% ...", round(j/nrow(mutations), 2) * 100))
  }
}
mutations$ref3_cod = ref3_cod
mutations$mut3_cod = mut3_cod
mutations$aachange = aachange
mutations$ntchange = ntchange
mutations$codonsub = codonsub
mutations$impact = impact
mutations$pid = sapply(RefCDS, function(x) x$protein_id)[mutations$geneind]
if (any(!is.na(wrong_ref))) {
  if (mean(!is.na(wrong_ref)) < 0.1) {
    warning(sprintf("%0.0f (%0.2g%%) mutations have a wrong reference base (see the affected mutations in dndsout$wrongmuts). Please identify the causes and rerun dNdScv.", 
                    sum(!is.na(wrong_ref)), 100 * mean(!is.na(wrong_ref))))
  }
  else {
    stop(sprintf("%0.0f (%0.2g%%) mutations have a wrong reference base. Please confirm that you are not running data from a different assembly or species.", 
                 sum(!is.na(wrong_ref)), 100 * mean(!is.na(wrong_ref))))
  }
  wrong_refbase = mutations[!is.na(wrong_ref), 1:5]
  mutations = mutations[is.na(wrong_ref), ]
}

# Annotate Indels with their change of coding region --------------------------
if (any(nrow(indels))) {
  indels = cbind(indels, 
                 data.frame(ref_cod = ".", mut_cod = ".", ref3_cod = ".", 
                            mut3_cod = ".", aachange = ".", ntchange = ".", 
                            codonsub = ".", impact = "no-SNV",
                            pid = sapply(RefCDS, 
                                         function(x) x$protein_id)[indels$geneind]))
  ins = nchar(gsub("-", "", indels$ref)) < nchar(gsub("-", "", indels$mut))
  del = nchar(gsub("-", "", indels$ref)) > nchar(gsub("-", "", indels$mut))
  multisub = nchar(gsub("-", "", indels$ref)) == nchar(gsub("-", "", 
                                                            indels$mut))
  l = nchar(gsub("-", "", indels$ref)) - nchar(gsub("-","", indels$mut))
  indelstr = rep(NA, nrow(indels))
  for (j in 1:nrow(indels)) {
    geneind = indels$geneind[j]
    pos = indels$start[j]:indels$end[j]
    if (ins[j]) {
      pos = c(pos - 1, pos)
    }
    pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, 
                      RefCDS[[geneind]]$strand)
    if (length(pos_ind) > 0) {
      inframe = (length(pos_ind)%%3) == 0
      if (ins[j]) {
        indelstr[j] = sprintf("%0.0f-%0.0f-ins%s", min(pos_ind), max(pos_ind),
                              c("frshift", "inframe")[inframe + 1])
      }
      else if (del[j]) {
        indelstr[j] = sprintf("%0.0f-%0.0f-del%s", min(pos_ind), max(pos_ind), 
                              c("frshift", "inframe")[inframe + 1])
      }
      else {
        indelstr[j] = sprintf("%0.0f-%0.0f-mnv", min(pos_ind), max(pos_ind))
      }
    }
  }
  indels$ntchange = indelstr
  annot = rbind(mutations, indels)
} else {
  annot = mutations
}
annot = annot[order(annot$sampleID, annot$chr, annot$pos), ]

# Compute number of expected synonymous, missense, etc mutations for every gene ----
message("[3] Estimating global rates...")
# to learn about definition of L matrix for each gene, see "SUPER IMPORTANT 
# NOTE" above
Lall = array(sapply(RefCDS, function(x) x$L), dim = c(192, 4, length(RefCDS)))
# The N matrix is very similar to the L matrix. It is created in the section
# "Annotate SNVs with their change of coding region". While L matrix tells you
# numbers of ALL possible mutations with a certain impact (i.e. missense) in
# the certain trinucleotide context, N tells you have many of those mutations
# actually happended in the gene in our cohort.
Nall = array(sapply(RefCDS, function(x) x$N), dim = c(192, 4, length(RefCDS)))
L = apply(Lall, c(1, 2), sum)
N = apply(Nall, c(1, 2), sum)

# the expected number of mutations of each type will be derived based on 
# poisson model (see function fit_substmodel above)
poissout = fit_substmodel(N, L, substmodel)
par = poissout$par
poissmodel = poissout$model
parmle = setNames(par[, 2], par[, 1])
mle_submodel = par
rownames(mle_submodel) = NULL
s1 = gsub("wmis", "wall", 
          gsub("wnon", "wall", gsub("wspl", "wall", substmodel)))
par1 = fit_substmodel(N, L, s1)$par
s2 = gsub("wnon", "wtru", gsub("wspl", "wtru", substmodel))
par2 = fit_substmodel(N, L, s2)$par
globaldnds = rbind(par, par1, par2)[c("wmis", "wnon", "wspl", "wtru", "wall"),]
sel_loc = sel_cv = NULL
genemuts = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), 
                      n_syn = NA, n_mis = NA, n_non = NA, n_spl = NA, 
                      exp_syn = NA, exp_mis = NA, exp_non = NA, exp_spl = NA,
                      stringsAsFactors = F)
genemuts[, 2:5] = t(sapply(RefCDS, function(x) colSums(x$N)))
mutrates = sapply(substmodel[, 1],
                  function(x) prod(parmle[base::strsplit(x, split = "\\*")[[1]]]))
genemuts[, 6:9] = t(sapply(RefCDS, function(x) colSums(x$L * mutrates)))
numrates = length(mutrates)
# if you look at genemuts table now, it will already have columns you're familiar with 
# gene_name n_syn n_mis n_non n_spl    exp_syn    exp_mis     exp_non     exp_spl

# Run dNdScv model without using covariates -----------------------------------
if (outp > 1) {
  message("[4] Running dNdSloc...")
  # see selfun_loc function above to see actual model
  sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), selfun_loc)))
  colnames(sel_loc) = c("wmis_loc", "wnon_loc", "wspl_loc", "pmis_loc",
                        "pall_loc")
  sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method = "BH")
  sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method = "BH")
  sel_loc = cbind(genemuts[, 1:5], sel_loc)
  sel_loc = sel_loc[order(sel_loc$pall_loc, sel_loc$pmis_loc, 
                          -sel_loc$wmis_loc), ]
}
nbreg = nbregind = NULL

# Run dNdScv model with using covariates --------------------------------------
# this section is very badly written, but essentially it fits the model with 
# using covariates and 
if (outp > 2) {
  message("[5] Running dNdScv...")
  if (is.null(cv)) {
    nbrdf = genemuts[, c("n_syn", "exp_syn")]
    model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf)
    message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", 
                    model$theta))
  }
  else {
    covs = as.matrix(covs[genemuts$gene_name, ])
    if (ncol(covs) > maxcovs) {
      warning(sprintf("More than %s input covariates. Only the first %s will be considered.", 
                      maxcovs, maxcovs))
      covs = covs[, 1:maxcovs]
    }
    nbrdf = cbind(genemuts[, c("n_syn", "exp_syn")], covs)
    if (nrow(genemuts) < mingenecovs) {
      model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf)
    } else {
      model = tryCatch({
        MASS::glm.nb(n_syn ~ offset(log(exp_syn)) +  ., data = nbrdf)
      }, warning = function(w) {
        MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf)
      }, error = function(e) {
        MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf)
      })
    }
    message(sprintf("    Regression model for substitutions (theta = %0.3g).", 
                    model$theta))
  }
  if (all(model$y == genemuts$n_syn)) {
    genemuts$exp_syn_cv = model$fitted.values
  }
  theta = model$theta
  nbreg = model

  sel_cv = as.data.frame(t(sapply(1:nrow(genemuts), selfun_cv)))
  colnames(sel_cv) = c("wmis_cv", "wnon_cv", "wspl_cv", 
                       "pmis_cv", "ptrunc_cv", "pallsubs_cv")
  sel_cv$qmis_cv = p.adjust(sel_cv$pmis_cv, method = "BH")
  sel_cv$qtrunc_cv = p.adjust(sel_cv$ptrunc_cv, method = "BH")
  sel_cv$qallsubs_cv = p.adjust(sel_cv$pallsubs_cv, method = "BH")
  sel_cv = cbind(genemuts[, 1:5], sel_cv)
  sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, 
                        sel_cv$ptrunc_cv, -sel_cv$wmis_cv), ]
  geneindels <- NULL
  if (nrow(indels) >= min_indels) {
    geneindels = as.data.frame(array(0, dim = c(length(RefCDS), 
                                                8)))
    colnames(geneindels) = c("gene_name", "n_ind", "n_induniq", "n_indused", 
                             "cds_length", "excl", "exp_unif", "exp_indcv")
    geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
    geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[, 1]])
    geneindels[is.na(geneindels[, 2]), 2] = 0
    geneindels$n_induniq = as.numeric(table(unique(indels[, -1])$gene)[geneindels[, 1]])
    geneindels[is.na(geneindels[, 3]), 3] = 0
    geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
    if (use_indel_sites) {
      geneindels$n_indused = geneindels[, 3]
    }
    else {
      geneindels$n_indused = geneindels[, 2]
    }
    geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
    min_bkg_genes = 50
    if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, 
                                                               "n_indused"]) == 0) {
      newkc = as.vector(sel_cv$gene_name[sel_cv$qallsubs_cv < 0.01])
      geneindels$excl = (geneindels[, 1] %in% newkc)
      if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, 
                                                                 "n_indused"]) == 0) {
        geneindels$excl = F
        message("    No gene was excluded from the background indel model.")
      }
      else {
        warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", 
                        paste(newkc, collapse = ", ")))
      }
    }
    geneindels$exp_unif = sum(geneindels[!geneindels$excl, 
                                         "n_indused"])/sum(geneindels[!geneindels$excl, 
                                                                      "cds_length"]) * geneindels$cds_length
    if (is.null(cv)) {
      nbrdf = geneindels[, c("n_indused", "exp_unif")][!geneindels[, 
                                                                   6], ]
      model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 
                             1, data = nbrdf)
      nbrdf_all = geneindels[, c("n_indused", "exp_unif")]
    }
    else {
      nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], 
                    covs)[!geneindels[, 6], ]
      if (sum(!geneindels$excl) < 500) {
        model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 
                               1, data = nbrdf)
      }
      else {
        model = tryCatch({
          MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + 
                         ., data = nbrdf)
        }, warning = function(w) {
          MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 
                         1, data = nbrdf)
        }, error = function(e) {
          MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 
                         1, data = nbrdf)
        })
      }
      nbrdf_all = cbind(geneindels[, c("n_indused", 
                                       "exp_unif")], covs)
    }
    message(sprintf("    Regression model for indels (theta = %0.3g)", 
                    model$theta))
    theta_indels = model$theta
    nbregind = model
    geneindels$exp_indcv = exp(predict(model, nbrdf_all))
    geneindels$wind = geneindels$n_indused/geneindels$exp_indcv
    geneindels$theta <- theta_indels
    geneindels$pind = pnbinom(q = geneindels$n_indused - 
                                1, mu = geneindels$exp_indcv, size = theta_indels, 
                              lower.tail = F)
    geneindels$qind = p.adjust(geneindels$pind, method = "BH")
    sel_cv = merge(sel_cv, geneindels, by = "gene_name")[, 
                                                         c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", 
                                                           "n_indused", "wmis_cv", "wnon_cv", "wspl_cv", 
                                                           "wind", "pmis_cv", "ptrunc_cv", "pallsubs_cv", 
                                                           "pind", "qmis_cv", "qtrunc_cv", "qallsubs_cv")]
    colnames(sel_cv) = c("gene_name", "n_syn", "n_mis", 
                         "n_non", "n_spl", "n_ind", "wmis_cv", "wnon_cv", 
                         "wspl_cv", "wind_cv", "pmis_cv", "ptrunc_cv", 
                         "pallsubs_cv", "pind_cv", "qmis_cv", "qtrunc_cv", 
                         "qallsubs_cv")
    sel_cv$pglobal_cv = 1 - pchisq(-2 * (log(sel_cv$pallsubs_cv) + 
                                           log(sel_cv$pind_cv)), df = 4)
    sel_cv$qglobal_cv = p.adjust(sel_cv$pglobal, method = "BH")
    sel_cv = sel_cv[order(sel_cv$pglobal_cv, sel_cv$pallsubs_cv, 
                          sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv), 
    ]
  }
}
if (!any(!is.na(wrong_ref))) {
  wrong_refbase = NULL
}
annot = annot[, setdiff(colnames(annot), c("start", "end", "geneind"))]
if (outmats) {
  dndscvout = list(globaldnds = globaldnds, sel_cv = sel_cv, 
                   sel_loc = sel_loc, annotmuts = annot, genemuts = genemuts, 
                   geneindels = geneindels, mle_submodel = mle_submodel, 
                   exclsamples = exclsamples, exclmuts = exclmuts, 
                   nbreg = nbreg, nbregind = nbregind, poissmodel = poissmodel,
                   wrongmuts = wrong_refbase, N = Nall, L = Lall)
} else {
  dndscvout = list(globaldnds = globaldnds, sel_cv = sel_cv, 
                   sel_loc = sel_loc, annotmuts = annot, genemuts = genemuts, 
                   geneindels = geneindels, mle_submodel = mle_submodel, 
                   exclsamples = exclsamples, exclmuts = exclmuts, 
                   nbreg = nbreg, nbregind = nbregind, poissmodel = poissmodel,
                   wrongmuts = wrong_refbase)
}