library("ArrayExpress")

GEOaccession = "GSE15765"
dataPath = file.path(getwd(), "data", GEOaccession)
supportFilesPrefix = "E-GEOD-15765"

mexp886 = list()
mexp886$path = dataPath
mexp886$rawFiles = list.files(dataPath, pattern="\\.CEL$")
mexp886$rawArchive = paste(supportFilesPrefix, ".raw.1.zip", sep="")
mexp886$processedFiles = paste(supportFilesPrefix, "-processed-data-1343526943.txt", sep="")
mexp886$processedArchive = paste(supportFilesPrefix, ".processed.1.zip", sep="")
mexp886$sdrf = paste(supportFilesPrefix, ".sdrf.txt", sep="")
mexp886$idf = paste(supportFilesPrefix, ".idf.txt", sep="")
mexp886$adf = "A-AFFY-37.adf.txt"

AEset = ae2bioc(mageFiles = mexp886)

colnames(pData(AEset))
fac = colnames(pData(AEset))[grep("Factor", colnames(pData(AEset)))]
fac

library("affyQCReport")
library("simpleaffy")

saqc=qc(AEset)
plot(saqc)

library("affyPLM")
dataPLM = fitPLM(AEset)

boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.22), outline = FALSE, col="lightblue", las=3, whisklty=0, staplelty=0)

# изменил масштаб, был -0.4; 0.4
Mbox(dataPLM, main="RLE", ylim = c(-0.8, 0.8), outline = FALSE, col="mistyrose", las=3, whisklty=0, staplelty=0)

library("arrayQualityMetrics")
arrayQualityMetrics(expressionset = AEset, outdir = "QAraw", force = TRUE, do.logtransform = TRUE, intgroup = fac)

