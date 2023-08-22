library(fst)
library(data.table)
library(Biostrings)

# Inputs ----------------------------------------------------------------------
mut_tab_path <- '/Users/maria/Desktop/BitBucket/LCCE_workshop_2023/figurecode/data/20221109_TRACERx421_mutation_table.fst'
signature_weigths_path <- '/Users/maria/Desktop/BitBucket/LCCE_workshop_2023/figurecode/data/20221110_TRACERx421_mutationSignature_weights.rds'
signature_profiles_path <- 'SBS_GRCh37.txt'

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
              family = poisson(link = log), weights = )
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
# Read TRACERx somatic driver mutations & signatures --------------------------
mutations <- read_fst(mut_tab_path)

# simplify mutation_id, please do not skip this step
muts_to_simplify <- grepl(';', mutations$mutation_id)
mutations[muts_to_simplify, ]$mutation_id = sapply(mutations[muts_to_simplify, ]$mutation_id,
                                                   function(x) sub(':', paste0('_', gsub('.*;', '', x), ':'), x))
mutations[muts_to_simplify, ]$mutation_id = gsub(';.*', '',
                                                 mutations[muts_to_simplify, ]$mutation_id)

# modify exonic.func column for clarity, please do not skip this step
mutations[grep('UNKNOWN', mutations$exonic.func), ]$exonic.func <- 'unknown'

# re-order columns for clarity, please do no skip this step
mutations <- mutations[, c('tumour_id', 'chr', 'start', 'ref', 'var')]
colnames(mutations) <- c('sampleID', 'chr', 'pos', 'ref', 'mut')
mutations$pos <- as.integer(mutations$pos)
mutations$chr <- gsub('chr', '', mutations$chr)

# read signatures weights and profiles
signature_weights = readRDS(signature_weigths_path)
signature_weights <- signature_weights$signature_weights_perTumour
signature_weights <- as.data.table(signature_weights)
setnames(signature_weights, 'tumour_id', 'sampleID')

signature_profiles <- fread('SBS_GRCh37.txt', sep = '\t', header = T,
                            stringsAsFactors = F)

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

# [STEP 1] Convert signatures profiles to 192 format --------------------------
# Data table signature_profiles was read from the file SBS_GRCh37.txt and 
# contains standard profiles of all signatures. This data table has 96 rows 
# each corresponding to the different trinucleotide context. First of all,
# change the format of how trinucleotides are written in the Type column so 
# they have have the same format as used in dNdScv table. For example:
# A[C>A]A will become ACA>AAA
# A[C>A]T will become ACT>AAT
# Second, rename column Type to column codon_change. 
# Third, convert signature_profiles data table to the data table with 192 rows
# (96 rows for the trinucleotide on forward strand and 96 rows for the 
# trinucleotides on reverse strand). 

# s_1A: convert column Type to the dNdScv format. 
# Tip: used function substr and paste
# ... insert your code here ...
# Assuming SBS_GRCh37 is a data.table
signature_profiles <- as.data.table(signature_profiles)
signature_profiles <- rbind(signature_profiles, signature_profiles)  #192 ROWS
# Modify the Type column directly using substr and paste
signature_profiles[, Type := paste0(substr(Type, 1, 1),
                            substr(Type, 3, 3),
                            substr(Type, 7, 7),
                            '>',
                            substr(Type, 1, 1),
                            substr(Type, 5, 5),
                            substr(Type, 7, 7))]

# s_1B: rename column Type to column codon_change
colnames(signature_profiles)[1] <- "codon_change"

# s_1C: convert signature_profiles data table to the data table with 192 rows
# (96 rows for the trinucleotide on forward strand and 96 rows for the 
# trinucleotides on reverse strand).
# Tip: to perform reverse compliment operation on trinucleotides, create 2 
# columns named ref3_cod and mut3_cod which will contain reference and mutated
# codon respectively. Use DNAString and reverseComplement functions from 
# Biostrings package to perform reverse compliment operation.
# !! Tip: make an assumption that strands contribute equally to signature profile

signature_profiles$ref3_cod <- substr(signature_profiles$codon_change,1,3)
signature_profiles$mut3_cod <- substr(signature_profiles$codon_change,5,7)
signature_profiles$ref3_cod_strand <-sapply(signature_profiles$ref3_cod,
                                    DNAString)
signature_profiles$mut3_cod_strand <- sapply(signature_profiles$mut3_cod,
                                     DNAString)
signature_profiles[97:192,]$ref3_cod_strand <-sapply(signature_profiles[97:192,]$ref3_cod_strand,
                                               reverseComplement)
signature_profiles[97:192,]$mut3_cod_strand <- sapply(signature_profiles[97:192,]$mut3_cod_strand,
                                               reverseComplement)
signature_profiles$ref3_cod_strand <-sapply(signature_profiles$ref3_cod_strand,
                                             as.character)
signature_profiles$mut3_cod_strand <- sapply(signature_profiles$mut3_cod_strand,
                                              as.character)
signature_profiles$codon_change <- paste0(signature_profiles$ref3_cod_strand,
                                          ">", signature_profiles$mut3_cod_strand)
#all weights divided by 2:make an assumption that strands contribute equally to signature profile
numeric_cols <- names(signature_profiles)[sapply(signature_profiles, is.numeric)]
signature_profiles[, (numeric_cols) := lapply(.SD, function(x) x / 2), .SDcols = numeric_cols]

# put the data table which you created in this section under name 
# signature_profiles_192
signature_profiles_192 <- signature_profiles

# [STEP 2]: Calculate distribution of mutations over 3nucl. context expected to be produced by signatures ----
# On this step calculate distribution of mutations' proportions over 
# trinucleotide context expected to be produced by signatures by multiplying
# signature weights by signature profiles

# s_2A: select from signature_profiles_192 signatures which are active in at least 1 patient.

#signature_weights should remove CRUK0418 first
signature_weights <- signature_weights[sampleID != "CRUK0418"]

weight_column_names <- colnames(signature_weights)[-1]
selected_SBS <- weight_column_names[colSums(signature_weights[, -1]) > 0]
rm(weight_column_names)

signature_weights <- signature_weights[, c("sampleID", selected_SBS),
                                       with = F]
signature_profiles_192 <- signature_profiles_192[, c("codon_change", selected_SBS),
                                                 with = F]

# s_2B: perform matrix multiplications with operator  %*%
#sample_idx <- unique(signature_weights$sampleID)
# put the data table which you created in this section under name 
# per_patient_trinucl_weigths. It should have columns sampleID plus one column
# per trinucleotide change. It should have 193 columns.

get_per_patient_trinucl_weigths <- function(sample_idx, 
                                             signature_weights,
                                             signature_profiles) {
  one_tumor_sw <- signature_weights[sampleID == sample_idx]
  one_tumor_sw <- one_tumor_sw[, colnames(one_tumor_sw) != 'sampleID', with = F]
  one_tumor_sw <- unlist(one_tumor_sw)
  # sp = Signature profile
  one_tumor_sp <- signature_profiles_192[, c('codon_change', names(one_tumor_sw)),
                                            with = F]
  one_tumor_sp <- as.data.frame(one_tumor_sp)
  rownames(one_tumor_sp) <- one_tumor_sp$codon_change
  one_tumor_sp <- one_tumor_sp[, colnames(one_tumor_sp) != 'codon_change']
  one_tumor_sp <- as.matrix(one_tumor_sp)
  one_tumor_expected <- data.table(trinucleotide = rownames(one_tumor_sp),
                                   per_trinucleotide_change = as.vector(one_tumor_sp %*%
                                                                 one_tumor_sw))
}
per_patient_trinucl_weigths <- lapply(unique(signature_weights$sampleID), 
                          get_per_patient_trinucl_weigths,signature_weights, signature_profiles_192)
per_patient_trinucl_weigths <- do.call(cbind,per_patient_trinucl_weigths)
per_patient_trinucl_weigths <- as.data.table(t(per_patient_trinucl_weigths))
per_patient_trinucl_weigths <- unique(per_patient_trinucl_weigths)
colnames(per_patient_trinucl_weigths) <- as.character(per_patient_trinucl_weigths[1,])
per_patient_trinucl_weigths <- per_patient_trinucl_weigths[-1,]
per_patient_trinucl_weigths <- cbind(signature_weights$sampleID, per_patient_trinucl_weigths)
colnames(per_patient_trinucl_weigths)[1] <- "sampleID"

# s_2C: sort the table in a way that the first column in it is sampleID and
# the other columns are in the order given in the names of trinucsubsind 
# vector. trinucsubsind is already defined.
# Sort the data table based on the order in trinucsubsind
per_patient_trinucl_weigths <- per_patient_trinucl_weigths[, c("sampleID", names(trinucsubsind)), with = FALSE]



# [STEP 3] create patient - specific N matrix ---------------------------------
# Reminder: in the classic version of dNdScv N matrix has 192 rows corresponding
# to the all possible 192 trinucleotide changes and 4 columns. Each column 
# corresponds to a one of possible impacts of mutation on CDS: first column - 
# synonymous, second column - nonsense, third column - missense, forth - 
# Essential splice site mutations. In the classic version of dNdScv number of 
# all synonymous mutations ACROSS ALL PATIENTS is computed and put in the 
# corresponding cell on N matrix. You goal in this section is to create the N
# matrix individual to each patient. (N for each patient)

# s_3A: extract columns sampleID, ref3_cod, mut3_cod and impact from mutations
# table and put it to the table named N_by_sampleID.
N_by_sampleID <- mutations[,c("sampleID", "ref3_cod", "mut3_cod","impact")]
# s_3B: in the table N_by_sampleID create a new column named codon_change which
# will hold a codon change code in dNdScv format (i.e. ACT>AAT). Use columns
# ref3_cod and mut3_cod to do so
N_by_sampleID$codon_change <- paste0(N_by_sampleID$ref3_cod,">",N_by_sampleID$mut3_cod)
# s_3C: remove from N_by_sampleID mutation which impact is not one of 
# Synonymous, Missense, Nonsense or Essential_Splice
N_by_sampleID <- N_by_sampleID[grepl("Synonymous|Missense|Nonsense|Essential_Splice", 
                                     N_by_sampleID$impact), ]
# s_3D: create in N_by_sampleID column with name impact_code which will be 
# equal to 1 if it's Synonymous mutation, 2 if it's Missense mutation, 3 if
# it's Nonsense mutation and 4 if it's Essential_Splice mutation.
# Create a function to map impact to impact_code
get_impact_code <- function(impact) {
  if (grepl("Synonymous", impact)) {
    return(1)
  } else if (grepl("Missense", impact)) {
    return(2)
  } else if (grepl("Nonsense", impact)) {
    return(3)
  } else if (grepl("Essential_Splice", impact)) {
    return(4)
  }
}

N_by_sampleID$impact_code <- sapply(N_by_sampleID$impact, get_impact_code)

# s_3E: restrict N_by_sampleID to columns sampleID, codon_change, impact_code
# and count number of mutation per each combination of sampleID, codon_change
# and impact_code.
N_by_sampleID <- as.data.table(N_by_sampleID)
N_by_sampleID_count <- N_by_sampleID[,.N,by=c("sampleID","codon_change",
                                              "impact_code")]

# s_3F: reshape table N_by_sampleID to the following format:
#        sampleID      codon_change     1  2  3  4
#        CRUK0001      AAA>AGA          1 NA NA NA
#        CRUK0001      AAA>ATA          NA  2 NA NA
# replace NA with 0
#dcast()
N_by_sampleID_count <- dcast(N_by_sampleID_count, 
                             sampleID + codon_change ~ impact_code, 
                             value.var = "N", fill = 0)


# s_3G: split table into list of tables by sampleID
N_by_sampleID_list <- split(N_by_sampleID_count,by="sampleID")
# s_3H: sort each table in a way that the order of trinucleotide changes in
# them are the same as in the names of trinucsubsind vector. trinucsubsind is
# already defined.

N_by_sampleID_list <- lapply(N_by_sampleID_list, function(table) {
  table <- table[match(names(trinucsubsind), table$codon_change), ]
  #table$codon_change <- names(trinucsubsind)
  return(table)
})

# s_3J: remove from each table column sampleID
# s_3JK: convert each table to matrix. Your matrixes should look this this:
#        [,1]   [,2]  [,3]    [,4]
#        [1,]    0    0    0    0
#        [2,]    1    0    0    0
#        [3,]    0    2    0    0
#        [4,]    0    0    0    0
# not all rows are printed 

N_by_sampleID_matrix <- lapply(N_by_sampleID_list, function(table) {
  table <- table[, !c("sampleID","codon_change")]
  table[is.na(table)] <- 0
  table_matrix <- as.matrix(table)
  return(table_matrix)
})

# by the end of this section you should have a list named N_by_sampleID with
# 431 elements, where each element is a mtrix 192x4.

####################### same result above
# [STEP 4] Modification of fit_substmodel function ----------------------------
# modify the template function below so that it takes one more additional 
# argument which should be named signatureW. This argument will be a vector of
# length 192 and will contain distribution of mutations' proportions over 
# trinucleotide context expected to be produced by signatures (see step 2).
# Think, if this vector needed to be modified? I.e. if we should submit here
# the original vector or maybe modify it somehow before?
# Check out documentation to function glm to see how weights can be
# incorporated
signatureW_matrix <- as.matrix(per_patient_trinucl_weigths[,-1])
rownames(signatureW_matrix) <- per_patient_trinucl_weigths$sampleID
signatureW_matrix <- apply(signatureW_matrix, MARGIN = 1, function(row) {
  row <- as.numeric(row)
  row <- 1 / row
  sum_row <- sum(row, na.rm = TRUE)
  row <- row / sum_row
  return(row)
})
signatureW_matrix <- t(signatureW_matrix)
colnames(signatureW_matrix) <- colnames(per_patient_trinucl_weigths[,-1])
signatureW_matrix <- as.data.frame(signatureW_matrix)
signatureW_matrix$sampleID <- rownames(signatureW_matrix)

# the same weight attribute to four impacts
signatureW_matrix <- apply(signatureW_matrix, MARGIN = 2, function(col) {
  rep(col, each = 4)
})
signatureW_matrix <- as.data.table(signatureW_matrix)
setcolorder(signatureW_matrix, 'sampleID')
setkey(signatureW_matrix, 'sampleID')
signatureW_matrix <- split(signatureW_matrix,by="sampleID")
signatureW_matrix <- lapply(signatureW_matrix,  t)
signatureW_matrix <- lapply(signatureW_matrix, function(table) {
  table <- table[match(names(trinucsubsind), rownames(table)), ]
  table_matrix <- as.matrix(table)
  return(table_matrix)
})


fit_substmodel_with_signatures = function(N, L, substmodel,signatureW) {
  l = c(L)
  n = c(N)
  w = c(signatureW)
  r = c(substmodel) 
  n = n[l != 0]
  r = r[l != 0]
  w = w[l != 0]
  l = l[l != 0]
  params = unique(base::strsplit(x = paste(r, collapse = "*"), 
                                 split = "\\*")[[1]])
  indmat = as.data.frame(array(0, dim = c(length(r), length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) {
    indmat[j, base::strsplit(r[j], split = "\\*")[[1]]] = 1
  }
  model = glm(formula = n ~ offset(log(l)) + . - 1, data = indmat, 
              family = poisson(link = log),weights = w)
  mle = exp(coefficients(model))
  ci = exp(confint.default(model))
  par = data.frame(name = gsub("`", "", rownames(ci)), mle = mle[rownames(ci)],
                   cilow = ci[, 1], cihigh = ci[, 2])
  return(list(par = par, model = model))
}

# [STEP 5] Compute number of expected synonymous, missense, etc mutations for every gene ----
# this is just a standard dNdScv way to compute L matrix - nothing changed here
message("[3] Estimating global rates...")
Lall = array(sapply(RefCDS, function(x) x$L), dim = c(192, 4, length(RefCDS)))
L = apply(Lall, c(1, 2), sum)

# s_5A: apply your modified function fit_substmodel_with_signatures to every
# matrix in list N_by_sampleID and save results in a list named 
# poissout_by_patient
poissout_by_patient <- lapply(names(N_by_sampleID_matrix),function(x){
  N <- N_by_sampleID_matrix[[x]]
  signatureW <- as.numeric(signatureW_matrix[[x]])
  result <- fit_substmodel_with_signatures(N = N, L = L, 
                                           substmodel=substmodel, 
                                           signatureW = signatureW)
  result$patient_id <- x
  return(result)
})


# s_5B: your list poissout_by_patient is actually list of lists. Extract
# from each element of the list an element named "par". These are the
# parameters of the model. One entry of a list should look like this:
#              name          mle        cilow     cihigh
#                 t 9.869420e-16 0.000000e+00        Inf
# `AAA>ACA` AAA>ACA 1.836426e+01 0.000000e+00        Inf
# `AAA>AGA` AAA>AGA 1.202515e+09 0.000000e+00        Inf
# `AAA>ATA` AAA>ATA 2.442667e+09 0.000000e+00        Inf
# not all rows are printed
par_by_patient <- lapply(poissout_by_patient, function(patient_result) {
  par_table <- patient_result$par
  return(par_table)
})
# s_5C: for each entry of the list, remove rows for which cilow is 0 and 
# cihigh is Inf. Tip: use function is.finite
par_by_patient <- lapply(par_by_patient, function(par_table) {
  par_table$mle[par_table$cilow == 0 | is.infinite(par_table$cihigh)] <- 0
  return(par_table)
})

# s_5D: for each entry of the list select only columns name and mle. patients in rows and mle in cols
par_by_patient <- lapply(par_by_patient, function(par_table) {
  par_table <- par_table[, !(names(par_table) %in% c("cilow", "cihigh")), drop = FALSE]
  return(par_table)
})

# s_5E: every entry of a list is patient specific. Now, sum up values in mle
# column of matrixes across all patients taking into account trinucleotide 
# context. Result of this should be a named vector which looks like this(sum up patients, get a vector 192)
#      AAA>AGA      AAA>ATA      AAC>ACC      AAC>AGC      AAC>ATC      AAG>ACG      AAG>AGG      AAG>ATG      AAT>ACT 
#      9.435630e-01 3.443699e-01 4.113114e+00 3.251753e+00 2.046960e+00 1.187832e+00 1.719782e+00 3.044378e+00 1.645949e+00 
#      AAT>AGT      AAT>ATT      ACA>AAA      ACA>AGA      ACA>ATA      ACC>AAC      ACC>AGC      ACC>ATC      ACG>AAG 
#      6.265443e+00 1.106065e+00 4.793340e+00 2.745451e+00 6.386596e+00 6.418332e+00 3.874873e+00 8.181891e+00 1.030336e+01 
# (not all values are printed)
# name this vector parmle
library(dplyr)
par <- bind_rows(par_by_patient)
par <- par %>%
  group_by(name) %>%
  summarise(mle = sum(mle))
parmle <- par$mle
names(parmle) <- par$name

#### continue running dNdScv as usual ####
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
  } else {
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
    } else {
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
      } else {
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
    } else {
      nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], 
                    covs)[!geneindels[, 6], ]
      if (sum(!geneindels$excl) < 500) {
        model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 
                               1, data = nbrdf)
      } else {
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



