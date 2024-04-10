require(kinship2)
library("readxl")
LizardSeqShort<-read_excel("../LeopardLizardSeqPedShort.xlsx") # Pedigree based on Sequoia (I used 1.0Extra here) exported to excel and supplemented with pedigree of ungenotyped individuals known from the zoo
LizardSeqFound<-read_excel("../LeopardLizardSeqPedFounder.xlsx") # Pedigree of only the founders

ped.lizardseqshort<-pedigree(id=LizardSeqShort$id, dadid=LizardSeqShort$sire, momid=LizardSeqShort$dam, sex=LizardSeqShort$sex)
LizardSeqFound<-LizardSeqFound[which(LizardSeqFound$id %in% c(1,2,3,4,5,19,20,'D1')),] #D1 = dummy, necessary; the other IDs are my known founder
ped.founder<-pedigree(id=LizardSeqFound$id, dadid=LizardSeqFound$sire, momid=LizardSeqFound$dam, sex=LizardSeqFound$sex,
                      affected=LizardSeqFound$affected, status=LizardSeqFound$status)

#Pedigree of founders, check to kinship matrix for degree of relationship
plot.pedigree(ped.lizardseqshort)
plot.pedigree(ped.founder, col="grey", cex=2, symbolsize=2)

# Pedigree from sequoia based on genotyped individual supplemented with known pedigree provided by the zoo from ungenotyped individuals, Since kinship 2 requires both dam and sire, unknown is replaced with Dummy parent
LizardSeqShort.status<-pedigree(LizardSeqShort$id, LizardSeqShort$sire, LizardSeqShort$dam, LizardSeqShort$sex,
                                affected=LizardSeqShort$affected, status=LizardSeqShort$status)
lizardseqshort.pedigree<-plot.pedigree(LizardSeqShort.status,col=ifelse(LizardSeqShort$inbreeding==0.0625, "green",
                                                        ifelse(LizardSeqShort$inbreeding==0.1250,"gold",
                                                               ifelse(LizardSeqShort$inbreeding==0.2500,"chocolate",
                                                                      ifelse(LizardSeqShort$inbreeding==0.3125,"red","grey")))))