library(GenomicRanges)
library(rtracklayer)
library(data.table)
files = list.files(".", pattern = "bedGraph$", recursive = TRUE, full.names = TRUE)
bodyMap = read.csv("/media/bowen_work/warwick/IHEC/bodyMap_withNewTissue.csv", header = TRUE, stringsAsFactors = FALSE)
bodyMap = rbind(bodyMap,
                rep("Plasma", ncol(bodyMap)),
                c(rep("testis", ncol(bodyMap)-1), "Endoderm"), c(rep("thymus", ncol(bodyMap)-1), "Endoderm"), c(rep("thyroid", ncol(bodyMap)-1), "Endoderm"))
colMap = list('Mesoderm' = as.matrix(c(231,138,195)), 'Endoderm' = as.matrix(c(141,160,203)),
              'Ectoderm' = as.matrix(c(0,204,102)), 'stem' = as.matrix(c(205,92,92)), 'Plasma' = as.matrix(c(255,215,0)))
colMap = lapply(colMap, function(x){ rownames(x) = c("red", "green", "blue") return(x)}
)
sample_data = data.frame(file = files,
                         patient = c('E12VA','E13SA', 'E13VA', 'E15VA', 'E18SA', 'E18VA', 'E23SA',
                                     'E23VA', 'E36VA', 'E40SA', 'E41SA', 'E42VA', 'E48SA'),
                         tissue = c('adipose', 'adipose', 'adipose', 'adipose', 'adipose', 'adipose',
                                    'adipose', 'adipose', 'adipose', 'adipose', 'adipose', 'adipose', 'adipose'), Sub_tissue = c('VisceralAdipocyte', 'SubcutaneousAdipocyte',
                                                                                                                                 'VisceralAdipocyte', 'VisceralAdipocyte', 'SubcutaneousAdipocyte', 'VisceralAdipocyte', 'SubcutaneousAdipocyte', 'VisceralAdipocyte', 'VisceralAdipocyte', 'SubcutaneousAdipocyte', 'SubcutaneousAdipocyte', 'VisceralAdipocyte', 'SubcutaneousAdipocyte'),
                         Health_status = c('Obese', 'Lean', 'Lean', 'Lean', 'Lean', 'Lean', 'Lean', 'Lean', 'Obese', 'Obese', 'Obese', 'Lean', 'Obese'))
for(i in files){
  bGraph = fread(i, data.table = FALSE)
  bGraph = makeGRangesFromDataFrame(bGraph, seqnames.field = "V1", start.field = "V2", end.field =
                                      "V3", keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE) #print("remove lambda, correct mcols and sort")
  tmp_score = mcols(bGraph)$V4
  mcols(bGraph) = NULL
  bGraph$score = tmp_score
  bGraph = bGraph[seqnames(bGraph) != 'lambda'] bGraph = sort(bGraph)
  tmp = sample_data[which(sample_data[,1] == i),]
  tline = new("TrackLine")
  tline@name=paste(tmp$tissue, tmp$Sub_tissue, tmp$Health_status, tmp$patient, sep = "_") tline@description=paste(tmp$tissue, tmp$Sub_tissue, sep = "_")
  col = colMap[[intersect(names(colMap), unique(bodyMap[grep(gsub("Derived","",tmp$tissue),
                                                             bodyMap$simple_tissue, ignore.case = TRUE),
                                                        "germ_layer"]))]]
  tline@color=as.integer(col)
  export.bedGraph(bGraph, "tmp.bedGraph", trackLine=tline) #print("cut out the extra useless bedgraph fields") system("cut -f 1-3,5 tmp.bedGraph > tmp_cut.bedGraph")
  outfile = paste("/media/scratch/lindsay/episcope/", tline@name, ".tdf", sep = '') cmd = paste("igvtools toTDF tmp_cut.bedGraph", outfile, "hg19")
  system(cmd)
  file.remove("tmp.bedGraph")
  file.remove("tmp_cut.bedGraph")
}