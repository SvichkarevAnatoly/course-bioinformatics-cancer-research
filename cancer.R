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
# dataPLM = fitPLM(AEset)
dataPLM = fitProbeLevelModel(AEset)

library("oligo")
oligo::NUSE(dataPLM) # посмотреть мануал для остальных параметров
oligo::RLE(dataPLM)

cAEset = oligo::rma(AEset)

groups = pData(AEset)[, fac] # убрать шум
groups[groups == "hepatocellular carcinoma"] = "HC"
groups[groups == "combined hepatocellular carcinoma and cholangiocarcinoma"] = "CHCAC"
groups[groups == "cholangiocarcinoma"] = "C"
groups

f = factor(groups)
f

design = model.matrix(~ 0 + f)
colnames(design)=sub("f", "", colnames(design))
design

library(limma)
conmat <- makeContrasts(HC-CHCAC, HC-C, CHCAC-C, levels = design)
colnames(conmat) <- c("HCvsCHCAC", "HCvsC", "CHCACvsC")
conmat

fit = lmFit(cAEset, design)
fitc <- contrasts.fit(fit, conmat)
fit2 = eBayes(fitc)

volcanoplot(fit2, coef="HCvsCHCAC", highlight=15)
volcanoplot(fit2, coef="HCvsC",     highlight=15)
volcanoplot(fit2, coef="CHCACvsC", highlight=15)

