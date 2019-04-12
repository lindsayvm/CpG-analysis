library(GenomicRanges) library(MethylSeekR) library(rtracklayer) library(data.table) library("BSgenome.Hsapiens.UCSC.hg19") sLengths = seqlengths(Hsapiens)
#load CpG islands
CpGislands.gr = import.bed("/media/bowen_work/warwick/Annotations/cpg_islands_hg19.bed") CpGislands.gr = suppressWarnings(resize(CpGislands.gr, 5000, fix = "center")) genome(CpGislands.gr) = NA
fn = list.files("/media/scratch/lindsay/episcope/", pattern = "*.bedGraph", full.names = TRUE) sapply(fn, function(X){
  meth.df = fread(X, data.table = FALSE)
  meth.gr = makeGRangesFromDataFrame(meth.df, seqnames.field = "V1", start.field = "V2", end.field = "V3", starts.in.df.are.0based = TRUE)
  meth.gr$T = rowSums(meth.df[ ,5:6]) #coverage meth.gr$M = meth.df[ ,5] #count
  ###UMRandLMR###
  #calculate FDR
  stats = calculateFDRs(m = meth.gr, CGIs = CpGislands.gr, num.cores = 40)
  FDR.cutoff = 5 #5%
  m.sel = 0.5
  n.sel = as.integer(names(stats$FDRs[as.character(m.sel), ][stats$FDRs[as.character(m.sel), ] <
                                                               FDR.cutoff])[1])
  n.sel
  #identify UMRs and LMRs
  meth.gr = meth.gr[seqnames(meth.gr) != "lambda"]
  UMRLMRsegments.gr = segmentUMRsLMRs(m = meth.gr, meth.cutoff = m.sel, nCpG.cutoff = n.sel, num.cores
                                      = 1, myGenomeSeq = Hsapiens, seqLengths = sLengths) #plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel) df = as.data.frame(UMRLMRsegments.gr)
  #remove all columns that...
  #colnames(df)
  df[4:11] = NULL
  write.table(df,file = gsub(".bedGraph$","_UMRsLMRs.bed", X), sep = "\t")
})