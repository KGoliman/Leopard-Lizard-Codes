# install.packages("kinship2") #only need to do once

# below are examples provided by the kinship tutorial (https://cran.r-project.org/web/packages/kinship2/kinship2.pdf)
# data(sample.ped)
# sample.ped
# sample.ped[1:10,]
# pedAll <- pedigree(id=sample.ped$id, 
#                   dadid=sample.ped$father, momid=sample.ped$mother, 
#                   sex=sample.ped$sex, famid=sample.ped$ped)
# print(pedAll)
# ped1basic <- pedAll['1']
# ped2basic <- pedAll['2']
# plot(ped1basic)
# plot(ped2basic)

require(kinship2) 
library("readxl") # not necessary, could download excel file as .csv and use read.csv
Lizard<-read_excel("../LeopardLizardR.xlsx") #provide information on ID, Dam, and Sire; "status" and "affected" are optional
Lizard<-Lizard[-which(Lizard$id %in% 121),] # 121 removed as it's newly introduced individual
Lizard<-Lizard[-which(Lizard$id %in% 124),] # ditto with 121, they had no relation to established individuals

ped.lizard<-pedigree(id=Lizard$id, dadid=Lizard$father, momid=Lizard$mother,
                     sex=Lizard$`Sex Type`)
plot.pedigree(ped.lizard)

# Lizard.alive<-subset(Lizard, Status=='Alive') # if only want data subset for currently living lizards
# Lizard.alive<-pedigree(id=Lizard.alive$id, dadid=Lizard.alive$father, momid=Lizard.alive$mother,
#                      sex=Lizard.alive$`Sex Type`)
# plot.pedigree(Lizard.alive) # only showing lizards currently alive

plot.pedigree(ped.lizard, col=ifelse(Lizard$Status=='Alive',"green",ifelse(Lizard$Status=='Dead',"red","blue")))
# Living = green, dead = red, released = blue

# My data uses Dead, Alive, and Released in the same column. Make new column for each with binary argument(1 = yes, 0 = no)
Lizard$dead<-with(Lizard, ifelse(Status == 'Dead', 1, 0))
Lizard$alive<-with(Lizard, ifelse(Status == 'Alive', 1, 0))
Lizard$released<-with(Lizard, ifelse(Status == 'Released', 1, 0))

Lizard.status<-pedigree(Lizard$id, Lizard$father, Lizard$mother,
                        Lizard$`Sex Type`, affected=Lizard$alive, status=Lizard$dead)
plot.pedigree(Lizard.status, col=ifelse(Lizard$dead,"black",1)) #Dead=crossed, Alive=filled, Released=Blank

#Finally, use the conditions of Lizard.status and pair it with previously determined inbreeding values to get pedigree
lizard.pedigree<-plot.pedigree(Lizard.status,col=ifelse(Lizard$Inbreeding==0.0625, "green",
                                        ifelse(Lizard$Inbreeding==0.1250,"gold",
                                         ifelse(Lizard$Inbreeding==0.2500,"chocolate",
                                          ifelse(Lizard$Inbreeding==0.3125,"red","grey")))))
