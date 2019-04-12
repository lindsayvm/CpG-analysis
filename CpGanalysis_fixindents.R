setwd("/media/scratch/lindsay/episcope/") library(GenomicRanges) library(rtracklayer)
library(data.table)
library(minfi)
library(limma)
library(gplots)
library(DMRcate)
library(ff)
##################Quick Load of files built below################## cpg.gr = import.bed("/media/scratch/lindsay/episcope/data/cpg.bed") ffload("/media/scratch/lindsay/episcope/data/meth_matrix") rownames(meth_ffdf) = mcols(cpg.gr)$name
meth_ffdf = meth_ffdf[ ,grep("KEP|precursor|AdiposeTissue_Unknown|CellLine|iPSC|ESC|Mesen|Celline|Hematopo|progenitor|early|[ g|G]lio|ancer|arcinoma|denoma|NOS", colnames(meth_ffdf), invert = TRUE)] load("/media/scratch/lindsay/episcope/data/cov.Rdata")
##################Add/remove files to/from IHEC##################
#all CpGs
cpg.df = fread("/media/bowen_work/warwick/Annotations/CpG_all_hg19.bed", data.table = FALSE) colnames(cpg.df) = c("seqnames", "start", "end", "name", "score", "extra")
cpg.gr = makeGRangesFromDataFrame(cpg.df, seqnames.field = "seqnames", start.field = "start", end.field = "end", starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE) export.bed(cpg.gr, con = "/media/scratch/lindsay/episcope/data/cpg.bed")
#Coverage matrix
fn = list.files("/media/scratch/lindsay/episcope/data/patients", pattern = "*.bedGraph", full.names = TRUE)
cov = rep(NA, length(cpg.gr))
for (i in fn){
  meth.df = fread(i, data.table = FALSE)
  meth.gr = makeGRangesFromDataFrame(meth.df, seqnames.field = "V1", start.field = "V2", end.field = "V3", starts.in.df.are.0based = TRUE)
  meth.gr$T = rowSums(meth.df[ ,5:6]) #coverage
  ov = as.matrix(findOverlaps(meth.gr, cpg.gr))
  #add scores of overlapping VA/SA ranges to new cov matrix,and NA for missing scores tmp_meth = rep(NA, length(cpg.gr))
  tmp_meth[ov[ ,2]] = mcols(meth.gr[ov[ ,1]])$T
  cov = cbind(cov, tmp_meth)
}
cov = cov[ ,2:ncol(cov)]
VASAnames = gsub(".*\\/|_.*", "", fn)
colnames(cov) = VASAnames
rownames(cov) = mcols(cpg.gr)$name
save(cov, file = "/media/scratch/lindsay/episcope/data/cov.Rdata")
##################Select hypo UMRs and hyper FMRs, exclude LMRs##################
#names of all samples bedfiles
VAfiles = list.files("/media/scratch/lindsay/episcope/data/patients/", pattern = "VA.*.UMRsLMRs.bed$", full.names = TRUE) #10
SAfiles = list.files("/media/scratch/lindsay/episcope/data/patients/", pattern = "SA.*.UMRsLMRs.bed$", full.names = TRUE) #3
otherfiles = list.files("/media/bowen_work/Data/IHEC/methylseekr/", pattern = "*.bed$", full.names = TRUE)
otherfiles = otherfiles[grep("AdiposeTissue_Unknown|CellLine|ESC", otherfiles, invert = TRUE)] #37
#List GRanges object of VA
VA_UMR = lapply(VAfiles,function(x){
  df = read.delim(x,stringsAsFactors = FALSE, header = TRUE, row.names)
  UMR.df = df[df$type == "UMR", ]
  UMR.gr = makeGRangesFromDataFrame(UMR.df, keep.extra.columns = TRUE)
  return(UMR.gr)
})
  names(VA_UMR) = gsub(pattern = ".*\\/|_.*", "", VAfiles ) #List GRanges object of other(SA)
  SA_UMR = lapply(SAfiles,function(x){
    df = read.delim(x,stringsAsFactors = FALSE, header = TRUE, row.names UMR.df = df[df$type == "UMR", ]
                    = 1)
    = 1)
  UMR.gr = makeGRangesFromDataFrame(UMR.df,
                                    return(UMR.gr)
  })
  names(SA_UMR) = gsub(pattern = ".*\\/|_.*", other_UMR = lapply(otherfiles, function(x){
    keep.extra.columns = TRUE)
    "", SAfiles )
  df = read.delim(x, stringsAsFactors = FALSE, header = TRUE) colnames(df) = c("seqnames","start","end","type")
  UMR.df = df[df$type == "UMR", 1:4]
  UMR.gr = makeGRangesFromDataFrame(UMR.df, keep.extra.columns = TRUE) return(UMR.gr)
  })
names(other_UMR) = gsub(pattern = ".*\\/|_.*", "", otherfiles ) other_UMR = c(SA_UMR, other_UMR)
#Merge all VA ranges in GR object
merge_VA_UMR = GRanges()
for(i in VA_UMR) merge_VA_UMR = c(i, merge_VA_UMR)
#Merge all VA and "other"(incl. SA) ranges in GR object
tmp = lapply(other_UMR, function(x){
  mcols(x)= NULL
  return(x) })
merge_other_UMR = GRanges()
for (i in tmp) merge_other_UMR = c(i, merge_other_UMR) tmp = merge_VA_UMR
mcols(tmp) = NULL
totalmerge_UMR = c(tmp, merge_other_UMR)
#disjoin VA,other(incl.SA) in non-overlapping regions totaldisjoin_UMR = disjoin(totalmerge_UMR)
#Ranges in totaldisjoin that NOT overlap with all "other"(incl.SA), thus hypo_overlap_UMR = subsetByOverlaps(totaldisjoin_UMR, merge_VA_UMR) hypo_VA_UMR = hypo_overlap_UMR[!hypo_overlap_UMR %over% merge_other_UMR] export.bed(hypo_VA_UMR, con = "hypo_VA_UMR_disjoin.bed")
#Ranges in totaldisjoin that overlap with all "others"(incl.SA) hyper_overlap_UMR = subsetByOverlaps(totaldisjoin_UMR, merge_other_UMR) hyper_VA_UMR = hyper_overlap_UMR[!hyper_overlap_UMR %over% merge_VA_UMR] export.bed(hyper_VA_UMR, con = "hyper_VA_UMR_disjoin.bed")
##################VA~SA UMR & LMR, for DMRcate################## #List GRanges object of VA
VA = lapply(VAfiles, function(x){
  are VA_unique
  df = read.delim(x, stringsAsFactors = FALSE, header = TRUE, row.names = 1) gr = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  return(gr)
})
names(VA) = gsub(pattern = ".*\\/|_.*", "", VAfiles ) #Merge all VA ranges in GR object
merge_VA = GRanges()
for(i in VA) merge_VA = c(i, merge_VA)
#List GRanges object of other(SA) SA = lapply(SAfiles, function(x){
SA = lapply(SAfiles, function(x){
  df = read.delim(x, stringsAsFactors = FALSE, header = TRUE, row.names = 1) gr = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  return(gr)
})
names(SA) = gsub(pattern = ".*\\/|_.*", "", SAfiles )
tmp = lapply(SA, function(x){
  mcols(x)= NULL
  return(x)
})
merge_SA = GRanges()
for (i in tmp) merge_SA = c(i, merge_SA) tmp = merge_VA
mcols(tmp) = NULL
Fat_totalmerge = c(tmp, merge_SA) #disjoin VA,SA in non-overlapping regions Fat_totaldisjoin = disjoin(Fat_totalmerge)
#Select all disjoined VA ranges from totaldisjoin
Fat_hypo_overlap = subsetByOverlaps(Fat_totaldisjoin, merge_VA) Fat_hypo_VA = Fat_hypo_overlap[!Fat_hypo_overlap %over% merge_SA] Fat_hypo_VA = reduce(Fat_hypo_VA)
export.bed(Fat_hypo_VA, con = "Fat_hypo_VA.bed")
#Ranges in totaldisjoin that overlap with all SA
Fat_hyper_overlap = subsetByOverlaps(Fat_totaldisjoin, merge_SA) Fat_hyper_VA = Fat_hyper_overlap[!Fat_hyper_overlap %over% merge_VA] Fat_hyper_VA = reduce(Fat_hyper_VA)
export.bed(Fat_hyper_VA, con = "Fat_hyper_VA.bed")

##################meth matrix overlap with UMR and FMRs##################
hypo_VA_UMR = import.bed("/media/scratch/lindsay/episcope/results/methylseekr/hypo_VA_UMR_disjoin.bed") 
hyper_VA_UMR = import.bed("/media/scratch/lindsay/episcope/results/methylseekr/hyper_VA_UMR_disjoin.bed")

#Reduce matrix to CpG of interest
VAhyperhypo_UMR.ls = list(hypo_VA_UMR, hyper_VA_UMR)
VAhyperhypo_UMR.gr = GRanges()
for(i in VAhyperhypo_UMR.ls) VAhyperhypo_UMR.gr = c(i, VAhyperhypo_UMR.gr)
VAhyperhypo_UMR.gr = sort(VAhyperhypo_UMR.gr)
export.bed(VAhyperhypo_UMR.gr, con = "/media/scratch/lindsay/episcope/results/methylseekr/VA_UMR_disjoin.bed") #disjoin
ov_UMR = as.matrix(findOverlaps(cpg.gr, VAhyperhypo_UMR.gr))
IHEC_meth_UMR = meth_ffdf[ov_UMR[ ,1], ]
#remove all rows and cols that are more than 80% NAs
IHEC_meth_UMR = IHEC_meth_UMR[which(!rowSums(is.na(IHEC_meth_UMR)) >= 0.8 * ncol(IHEC_meth_UMR)), which(!colSums(is.na(IHEC_meth_UMR)) >= 0.8 * nrow(IHEC_meth_UMR))]
dim(IHEC_meth_UMR)
IHEC_meth_UMR = as.matrix(IHEC_meth_UMR)
save(IHEC_meth_UMR, file = "/media/scratch/lindsay/episcope/results/intersect/IHEC_meth_UMR_ov.RData")


##################Limma##################
IHEC_meth_UMR = logit2(IHEC_meth_UMR / 100)
#Replace all inf with NA
IHEC_meth_UMR[which(!is.finite(IHEC_meth_UMR))] = NA
IHEC_meth_UMR = IHEC_meth_UMR[rowSums(is.na(IHEC_meth_UMR)) != ncol(IHEC_meth_UMR), colSums(is.na(IHEC_meth_UMR)) != nrow(IHEC_meth_UMR)]
dim(IHEC_meth_UMR)
save(IHEC_meth_UMR, file = "/media/scratch/lindsay/episcope/results/limma/IHEC_UMR_m.RData") #make design
VA_cols = colnames(IHEC_meth_UMR)[grep("Visceral|omentum", colnames(IHEC_meth_UMR))]
type = as.factor(ifelse(colnames(IHEC_meth_UMR) %in% VA_cols, "VA", "Other"))
design = model.matrix( ~ 0 + type)

cm = makeContrasts(c = typeVA-typeOther, levels = design)
#comparison
lfit = lmFit(IHEC_meth_UMR, design) #not robust, bcs VA is small dataset
cfit = contrasts.fit(lfit, cm)
efit = eBayes(cfit) #Bayes ttest
#Summary
tt_UMR = topTable(efit, coef = "c", number = nrow(IHEC_meth_UMR))
save(tt_UMR, file = "/media/scratch/lindsay/episcope/results/limma/IHEC_UMR_lm.RData") #Criteria: neg FC CpG within UMR
tt_hypo = tt_UMR[which(tt_UMR$logFC < 0), ]
dim(tt_hypo)
start = as.integer(gsub(".*_", "", rownames(tt_hypo)))
df = data.frame("seqnames" = gsub("_.*", "", rownames(tt_hypo)), start, "end" = start + 1) df$name = rownames(tt_hypo)
gr = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
hypo = subsetByOverlaps(gr, hypo_VA_UMR)
length(hypo)
#pos FC CpG within FMR
tt_hyper = tt_UMR[which(tt_UMR$logFC > 0), ]
dim(tt_hyper)
start = as.integer(gsub(".*_", "",rownames(tt_hyper)))
df = data.frame("seqnames" = gsub("_.*", "", rownames(tt_hyper)), start, "end" = start + 1) df$name = rownames(tt_hyper)
gr = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
hyper = subsetByOverlaps(gr, hyper_VA_UMR)
length(hyper)
VA_SAother = sort(c(hypo, hyper))
length(VA_SAother)
##################Difference in means##################
#Select sign in beta values
load(file = "/media/scratch/lindsay/episcope/results/intersect/IHEC_meth_UMR_ov.RData")#beta values tt_names = tt_UMR[which(rownames(tt_UMR) %in% mcols(VA_SAother)$name), ] #only neg CpGs in UMRs and pos in FMR
Pval = tt_names[which(tt_names$adj.P.Val <= 0.05), ] #for 0.01: 11014 #but after other filterting steps, no diff between 01 or 05

#sum(rownames(Pval) %in% mcols(hyper)$name)
IHEC_meth_UMR = IHEC_meth_UMR[rownames(IHEC_meth_UMR) %in% rownames(Pval), ]
dim(IHEC_meth_UMR)
#VAgroup:
VA_cols = colnames(IHEC_meth_UMR)[grep("Visceral|omen", colnames(IHEC_meth_UMR))]
VA_UMR = IHEC_meth_UMR[ ,colnames(IHEC_meth_UMR) %in% VA_cols]
VA_UMR_means = rowMeans(VA_UMR, na.rm = TRUE)
#othergroup
other_cols = colnames(IHEC_meth_UMR)[grep("Visceral|omentum", colnames(IHEC_meth_UMR), invert = TRUE)] other_UMR = IHEC_meth_UMR[ ,!colnames(IHEC_meth_UMR) %in% VA_cols]
other_UMR_means = rowMeans(other_UMR, na.rm = TRUE)
#absolute difference in means
UMR_means = abs(VA_UMR_means - other_UMR_means)
sum(UMR_means >= 50, na.rm = TRUE)
DE_UMR = IHEC_meth_UMR[which(UMR_means >= 50), ]
DE_UMR = DE_UMR[which(!rowSums(is.na(DE_UMR)) >= 0.8 * ncol(DE_UMR)), which(!colSums(is.na(DE_UMR)) >= 0.8 * nrow(DE_UMR))]
dim(DE_UMR)
save(DE_UMR, file = "/media/scratch/lindsay/episcope/results/means/DE_VA_UMR_50.RData")
##################DMRcate##################
#VA~SA URMLMR
#make design #merge_VA and merge_SA; add count and coverage columns
fn = list.files("/media/scratch/lindsay/episcope/data/patients", pattern = "*CpG.bedGraph", full.names = TRUE)
samples = sapply(fn, function(y){
  CpGs.df = fread(y, data.table = FALSE)
  colnames(CpGs.df) = c("chr","pos","V3", "V4", "V5", "V6") CpGs.df$N = rowSums(CpGs.df[,5:6]) #coverage, total count CpGs.df$X = CpGs.df[,5] #count; methylated count return(CpGs.df[,c("chr", "pos", "N", "X")])
}, simplify = FALSE)
obj_bsseq = makeBSseqData(samples, colnames(Fat_meth))
save(obj_bsseq, file = "bsseq.RData")
#call differentially methylated CpG sites
VA_cols = colnames(Fat_meth)[grep("VA|omentum",colnames(Fat_meth))]
SA_cols = colnames(Fat_meth)[grep("SA|subcut",colnames(Fat_meth))]
DSSres = DMLtest(obj_bsseq, group1 = VA_cols , group2 = SA_cols , smoothing = FALSE) save(DSSres, "DSSres.RData")
wgbsannot = cpg.annotate("sequencing", DSSres)
wgbs.DMRs = dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1) save(wgbs.DMRs, file = "wgbsDMR.RData")
wgbs.ranges= extractRanges(wgbs.DMRs, genome = "hg19")
save(wgbs.ranges, gsub("wgbsranges.RData"))
##################limma DMRcate overlap##################
#wgbsranges is DMRcate ouput from VA_SA UMR,LMR
#DE_VA_UMR_50 only contains UMRs and will therefore also exclude LMRs load("/media/scratch/lindsay/episcope/results/means/DE_VA_UMR_50.RData") start = as.integer(gsub(".*_", "", rownames(DE_UMR)))
df = data.frame("seqnames" = gsub("_.*", "", rownames(DE_UMR)), start, "end" = start + 1) gr = makeGRangesFromDataFrame(df)
VA_SAother = subsetByOverlaps(gr, wgbs.ranges)
# select count data from matrix for ranges of interest (single CpG sites)
ov_UMR = as.matrix(findOverlaps(cpg.gr, VA_SAother))
DE_UMR = meth_ffdf[ov_UMR[ ,1], ]
# NA check for only samples, this is already done for CpGs
DE_UMR = DE_UMR[ ,which(!colSums(is.na(DE_UMR)) >= 0.80 * nrow(DE_UMR))] # 741 254, of which 383 hypo and 358 hyper
COI = cpg.gr[ov_UMR[ ,1]]
VA_UMR = import.bed("/media/scratch/lindsay/episcope/results/methylseekr/VA_UMR_disjoin.bed") #disjoin
ROI = (subsetByOverlaps(VA_UMR, COI))
save(ROI, file = "ROI_limmaDMR.RData")
##################Sliding window##################
#data = numeric vector, window =how many CpG(window size only takes odd numbers),(step = how many CpGs to slide along)
slide = function(data, window, step){
  total = length(data)
  endPos = (window / 2) - 0.5
  spots = seq(from = (1 + endPos), to = (total - endPos), by = step) result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = mean(data[(spots[i] - endPos):(spots[i] + endPos)], na.rm = TRUE) }
  return(result)
}
####select from meth_ffdf all CpGs within ROI
#assign names to each CpG
names.df = as.data.frame(cpg.gr)
mcols(cpg.gr)$name = paste(names.df$seqnames, names.df$start, sep = "_") names.df = as.data.frame(ROI)
mcols(ROI)$name = paste(names.df$seqnames, names.df$start, sep = "_") #select all CpGs within ROI
tmp.gr = resize(ROI, 100 * width(ROI), fix = "center")
cpg_cutdown = subsetByOverlaps(cpg.gr, tmp.gr)
#select cpg_cutdown from met_ffdf with all CpGs  

methOI = meth_ffdf[mcols(cpg_cutdown)$name, ]
dim(methOI)
VA_cols = colnames(methOI)[grep("Visceral|omentum", colnames(methOI))]
other_cols = colnames(methOI)[grep("Visceral|omentum", colnames(methOI), invert = TRUE)] windowsize = 5
slideDiff_all = numeric()
for(i in 1:length(ROI)){
  print(i)
  ROI_tmp = ROI[i]
  COI_tmp = subsetByOverlaps(cpg_cutdown, ROI_tmp) #Get all CpGs that fall into region of interest COI_tmp = sort(COI_tmp) #genomic order
  COI_names = mcols(COI_tmp)$name
  j=1
  while(j < (windowsize / 2) + 0.5){
    #follow gives preceding nearest range neighbour and precede the following
    #adds CpG dummies to start and end of GR object, to calculate a near as possible mean follow = cpg_cutdown[follow(COI_tmp[1], cpg_cutdown)]
    precede = cpg_cutdown[precede(COI_tmp[length(COI_tmp)], cpg_cutdown)]
    COI_tmp = c(follow, COI_tmp, precede)
    COI_tmp = sort(COI_tmp)
    j=j+1
  }
  #only select the COI of methOI (COI is not resized)
  methOI_tmp = methOI[mcols(COI_tmp)$name, ]
  VAMethMean = rowMeans(methOI_tmp[ ,VA_cols], na.rm = TRUE) otherMethMean_all = rowMeans(methOI_tmp[ ,other_cols], na.rm = TRUE) VAslide = slide(VAMethMean, windowsize, 1)
  otherslide_all = slide(otherMethMean_all, windowsize, 1) slideDiff_all_tmp = abs(VAslide - otherslide_all) names(slideDiff_all_tmp) = COI_names
  slideDiff_all = c(slideDiff_all, slideDiff_all_tmp)
}
#save(slideDiff_all, file = "/media/scratch/lindsay/episcope/results/slidingwindow/sliding_disjoin_new.Rdata") hist(slideDiff_all, breaks = 100, main = "Difference in Mean of CpGs", xlab = "mean_VA - mean_other") #absolute difference in mean of CpGs
target_sites_UMR = names(slideDiff_all[which(slideDiff_all >= 50)])
target_sites_UMR = cpg_cutdown[which(mcols(cpg_cutdown)$name %in% target_sites_UMR)]
target_sites_UMR = subsetByOverlaps(VA_UMR, target_sites_UMR)
export.bed(target_sites_UMR, con = "/media/scratch/lindsay/episcope/results/slidingwindow/targetsites_UMR_disjoin_w5.bed") target_sites_UMR = import.bed("/media/scratch/lindsay/episcope/results/slidingwindow/targetsites_UMR_disjoin_w5.bed")
##################AmpliconCharacteristics################## ROI = target_sites_UMR
#Amplicon Feature filtering
#at least 80 nc long
# hist(width(ROI), breaks = 100) # range(width(ROI)) # table(width(ROI) > 80) ROI = ROI[width(ROI) > 80]
#at least 10 CpGs
tmp = subsetByOverlaps(cpg.gr, ROI)
c_count = sapply(ROI, function(x) sum(tmp %over% x)) #set cutoff on minimal 10 CpG # plot(c_count, width(ROI))
ROI = ROI[which(c_count > 10)]
# #more than 1% CpGs
# tmp = subsetByOverlaps(cpg.gr, ROI)
# c_count = sapply(ROI, function(x) sum(tmp %over% x))
# ROI = ROI[which(c_count > 0.01 * width(ROI))]
export.bed(ROI, con = "/media/scratch/lindsay/episcope/results/ROI_w5.bed")
#check for equal partition of hyper and hypo df_hyper = as.data.frame(hyper_VA_UMR)

df_hyper$name = paste(df_hyper$seqnames, df_hyper$start, sep = "_") df_hypo = as.data.frame(hypo_VA_UMR)
df_hypo$name = paste(df_hypo$seqnames, df_hypo$start, sep = "_") df_ROI = as.data.frame(ROI)
df_ROI$name = paste(df_ROI$seqnames, df_ROI$start, sep = "_") sum(df_ROI$name %in% df_hyper$name)
sum(df_ROI$name %in% df_hypo$name)
##################After being manually modified in excel################## ROI.df = read.csv("/media/bowen_work/Lindsay/ROI_w5_cutdown.csv")
ROI.gr = makeGRangesFromDataFrame(ROI.df, keep.extra.columns = TRUE) ROI.gr = sort(ROI.gr)
export.bed(ROI.gr, con = "/media/scratch/lindsay/episcope/results/ROI_w5_cutdown.bed")
#Create extended bed file
write.table(ROI.df, "/media/scratch/lindsay/episcope/ROI_w5_cutdown.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#second primer pool (Primersuite)
ROI.df = read.csv("/media/bowen_work/Lindsay/ROI_excel/pool2.csv") ROI.gr = makeGRangesFromDataFrame(ROI.df, keep.extra.columns = TRUE) ROI.gr = sort(ROI.gr)
##################Visualize##################
#Boxplots: each region on 1 pdf page with detail information
VA_cols = colnames(meth_ffdf)[grep("Visceral|omentum", colnames(meth_ffdf))]
other_cols = colnames(meth_ffdf)[!colnames(meth_ffdf) %in% VA_cols]
#other_cols = other_cols[grep("[S|s]ubc", other_cols, invert = TRUE)] #NO improvement
rm(list = setdiff(ls(), c("meth_ffdf", "cpg.gr", "VA_cols", "other_cols", "ROI" , "cov"))) pdf("/media/bowen_work/Lindsay/figures/boxplots/pool2.pdf", paper = "a4r", width = 11, height = 7) for(i in 1:length(ROI)){
  print(i)
  COI = subsetByOverlaps(cpg.gr, ROI[i]) #all cpgs wihtin ROI
  COI = sort(COI)
  #Coverage CpG cutoff (adipose)
  cov_COI = cov[COI$name, ]
  cov_COI = cov_COI[which(!rowSums(is.na(cov_COI)) >= 0.6 * ncol(cov_COI)), ] #assumption that cov is
  low when >60% is NA
  cov_COI = cov_COI[which(rowMeans(cov_COI, na.rm = TRUE) > 10), ] if(sum(COI$name %in% rownames(cov_COI)) != 0){
    COI = COI[COI$name %in% rownames(cov_COI)] #all cpgs wihtin ROI
    #methylation of CpGs wihtin ROI(high coverage) methOI_VA = meth_ffdf[COI$name, VA_cols] methOI_Other = meth_ffdf[COI$name, other_cols] boxplot(t(methOI_VA), at = 1:nrow(methOI_VA) *
    ylim = c(0, 110), start = 0, pch = 18, paste(as.character(i),":", COI$name[1], "-", COI$name[length(COI)]))
boxplot(t(methOI_Other), at = 1:nrow(methOI_VA) * 2 + 1, col = "blue", add = TRUE, pch = 18, xaxt = "n", yaxt = "n" )
legend("topleft", c("VA", "Other"), text.col = c("red", "blue")) }
}
dev.off() #return to null device
#Boxplots: all 75 regions on 2 pages for Appendix
VA_cols = colnames(meth_ffdf)[grep("Visceral|omentum", colnames(meth_ffdf))]
other_cols = colnames(meth_ffdf)[!colnames(meth_ffdf) %in% VA_cols]
#other_cols = other_cols[grep("[S|s]ubc", other_cols, invert = TRUE)] #NO improvement pdf("/media/bowen_work/Lindsay/figures/boxplots/pool2.pdf", paper = "a4r", width = 11, height = 7) for(i in 1:length(ROI)){
print(i)
COI = subsetByOverlaps(cpg.gr, ROI[i]) #all cpgs wihtin ROI COI = sort(COI)
#Coverage CpG cutoff (adipose)

cov_COI = cov[COI$name, ]
cov_COI = cov_COI[which(!rowSums(is.na(cov_COI)) >= 0.6 * ncol(cov_COI)), ] #assumption that cov is low when >60% is NA
cov_COI = cov_COI[which(rowMeans(cov_COI, na.rm = TRUE) > 10), ] if(sum(COI$name %in% rownames(cov_COI)) != 0){
  COI = COI[COI$name %in% rownames(cov_COI)] #all cpgs wihtin ROI
  #methylation of CpGs wihtin ROI(high coverage) methOI_VA = meth_ffdf[COI$name, VA_cols] methOI_Other = meth_ffdf[COI$name, other_cols] boxplot(t(methOI_VA), at = 1:nrow(methOI_VA) *
  ylim = c(0, 110), start = 0, pch = 18, paste(as.character(i),":", COI$name[1], "-", COI$name[length(COI)]))
boxplot(t(methOI_Other), at = 1:nrow(methOI_VA) * 2 + 1, col = "blue", add = TRUE, pch = 18, xaxt = "n", yaxt = "n" )
legend("topleft", c("VA", "Other"), text.col = c("red", "blue")) }
}
dev.off() #return to null device
##################Visualize##################
ov = as.matrix(findOverlaps(cpg.gr, ROI.gr))
DE = as.matrix(meth_ffdf[ov[ ,1], ])
DE = DE[which(!rowSums(is.na(DE)) >= 0.6 * ncol(DE)), which(!colSums(is.na(DE)) >= 0.6 * nrow(DE))]
VA_cols = colnames(DE)[grep("Visceral|omentum", colnames(DE))] SA_cols = colnames(DE)[grep("[s|S]ubcut", colnames(DE))] Muscle_cols = colnames(DE)[grep("Muscle", colnames(DE))] Breast_cols = colnames(DE)[grep("Breast", colnames(DE))] Liver_cols = colnames(DE)[grep("Liver", colnames(DE))]
#MDS
groups = rep("Other", ncol(DE)) groups[colnames(DE) %in% VA_cols] = "Visceral" groups[colnames(DE) %in% SA_cols] = "Subcutaneous" groups[colnames(DE) %in% Muscle_cols] = "Muscle" groups[colnames(DE) %in% Breast_cols] = "Breast" groups[colnames(DE) %in% Liver_cols] = "Liver" groups = as.factor(groups)
names = gsub("_.*", "", colnames(DE))
mdsPlot(DE, sampGroups = groups, sampNames = names,
        pal = c("purple", "yellow2", "dodgerblue", "grey40", "red2", "green2"),
        main = "VA ~ Reference Tissue", legendPos = "topleft", legendNCol = 2, pch = 17)
#Heatmap
cols = rep("black", ncol(DE))
cols[colnames(DE) %in% VA_cols] = "green2"
cols[colnames(DE) %in% SA_cols] = "red"
cols[colnames(DE) %in% Muscle_cols] = "dodgerblue"
cols[colnames(DE) %in% Breast_cols] = "purple2"
cols[colnames(DE) %in% Liver_cols] = "yellow2"
names = gsub("_.*", "", colnames(DE))
heatmap.2(DE, col = bluered, ColSideColors = cols, trace = "none", density.info = "none", na.color = "grey80",
          labCol = names) #, Rowv = FALSE) #dendrogram = "column"
##################miniSeq##################
fn = list.files("/media/scratch/lindsay/episcope/miniSeq/", pattern = "*.bedGraph", full.names = TRUE, recursive = TRUE)
adi.fn = fn[grep("sub|vis|adi", fn)]
#make ROI of the right 29 primers
ROI.df = read.csv("/media/bowen_work/Lindsay/primeramplicon/ROI_excel/ROIprimers28.csv")

ROI.gr = makeGRangesFromDataFrame(ROI.df, keep.extra.columns = TRUE)
ROI.gr = sort(ROI.gr)
COI.gr = subsetByOverlaps(cpg.gr, ROI.gr) #Get all CpGs that fall into region of interest COI.gr = sort(COI.gr) #genomic order
adi.mat = matrix(rep(NA, length(COI.gr))) for (i in adi.fn){
  new.df = fread(i, data.table = FALSE)
  new.gr = makeGRangesFromDataFrame(new.df, seqnames.field = "V1", start.field = "V2", end.field = "V3", starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  #matrix all cpgs
  ov = as.matrix(findOverlaps(new.gr, COI.gr)) tmp_meth = rep(NA, length(COI.gr))
  tmp_meth[ov[ ,2]] = new.df[ov[ ,1], "V4"] #V4=score adi.mat = cbind(adi.mat, tmp_meth)
}
adi.mat = adi.mat[ ,2:4]
adi.names = gsub(".*\\/|_.*", "", adi.fn) 
colnames(adi.mat) = adi.names 
rownames(adi.mat) = COI.gr$name
#meth.mat = as.matrix(meth.mat)
adi.mat = adi.mat[rowSums(is.na(adi.mat)) != ncol(adi.mat), colSums(is.na(adi.mat)) != nrow(adi.mat)] adi.vc = rowMeans(adi.mat, na.rm = TRUE)
VASA_cols = colnames(meth_ffdf)[grep("Visceral|omentum|[s|S]ubcut", colnames(meth_ffdf))] VASA.mat = meth_ffdf[rownames(adi.mat), VASA_cols]
VASA.vc = rowMeans(VASA.mat, na.rm = TRUE)
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(VASA.vc, meth.vc, colramp = jet.colors, nrpoints = 0, xlab = "WGBS", ylab = "miniSeq") abline(lm(meth.vc~VASA.vc), col = "black", lwd = 3, lty = "dotted")
text(1, 0, "y = 12 * 0.9x", col = "white", adj = c(-3.5,-5), lwd = 5)
##################perc barplots################## meth.fn = fn[grep("perc", fn)]
meth.mat = matrix(rep(NA, length(COI.gr)))
for (i in meth.fn){
  new.df = fread(i, data.table = FALSE)
  new.gr = makeGRangesFromDataFrame(new.df, seqnames.field = "V1", start.field = "V2", end.field = "V3", starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  #matrix all cpgs
  ov = as.matrix(findOverlaps(new.gr, COI.gr)) tmp_meth = rep(NA, length(COI.gr))
  tmp_meth[ov[ ,2]] = new.df[ov[ ,1], "V4"] #V4=score meth.mat = cbind(meth.mat, tmp_meth)
}
rownames(meth.mat) = COI.gr$name
meth.mat = meth.mat[rownames(adi.mat) ,2:4] meth.names = gsub(".*\\/|_.*", "", meth.fn) colnames(meth.mat) = meth.names
meth.vc = colMeans(meth.mat, na.rm = TRUE)
x = c(meth.vc["0perc"], meth.vc["50perc"], meth.vc["100perc"]) barplot(x, ylim = c(0,100))
##################450k data##################
HM450 = readRDS("/media/work/Data/Public/HM450/grset_2017-09-04.Rds") gr = granges(HM450)
df = as.data.frame(gr)
df$name = paste(df$seqnames, df$start, sep = "_")
gr = granges(df)


gr = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
pheno.df = read.csv("/media/bowen_work/Lindsay/HM450/phenodata_2017-09-04.csv", row.names = 1) adi.df = pheno.df[grep("adipose", pheno.df$tissue), grep("tissue|pheno|cancer_state", colnames(pheno.df))]
adi.df = adi.df[adi.df$cancer_state == "normal", ]
betas = getBeta(HM450) #hit twice for dependencies to load Mbetas = (as.matrix(betas)) * 100 #485512 5472
ov = as.matrix(findOverlaps(gr, ROI.gr)) #42 cpgs in 34 regions of pool2
DE = as.matrix(Mbetas[ov[ ,1], ])
DE = DE[which(!rowSums(is.na(DE)) >= 0.6 * ncol(DE)), which(!colSums(is.na(DE)) >= 0.6 * nrow(DE))] DE_adi = DE[ ,rownames(adi.df)]
DE_ref = DE[ ,grep("adipose", pheno.df$tissue, invert = TRUE)]
#DE_ref_samp = DE_ref[ ,sample(colnames(DE_ref), size = 1000)]
DE = cbind(DE_adi, DE_ref)


















  
  
  