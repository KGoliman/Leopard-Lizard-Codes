library(sequoia)
GenoFile<-GenoConvert(InFile="../inputSequoia_2.raw.raw", InFormat = "raw") #There are  95  individuals and  615  SNPs.
GenoStats<-SnpStats(GenoFile)
View(GenoStats) #Check for MAF and missingness qualities

LizardPed<-read.csv("../LeopardLizardPedShorter.csv", header = TRUE) #Pedigree information: ID, Sex, Dam, Sire
LizardLH<-read.csv("../LeopardLizardLHShorter.csv", header = TRUE) #Life-history information: ID, Sex, BirthYear, & estimated Min & Max BirthYear if unknown
Geno.checked<-CheckGeno(GenoFile, Return = "GenoM")

#Check for call rate, higher = better, if any are too low, e.g. <50%, may influence data, consider dropping low quality SNPs based on your standards
SnpCallRate <- apply(GenoFile, MARGIN=2,  FUN = function(x) sum(x!=-9)) / nrow(GenoFile)
IndivCallRate <- (apply(GenoFile, MARGIN = 1, FUN = function(x) sum(x != -9)) + rnorm(nrow(GenoFile), mean = 0, sd = 0.01)) / ncol(GenoFile)
hist(SnpCallRate, breaks=80, col="tan")
goodsamples <- rownames(GenoM)[SnpCallRate > 0.95]
badsamples <- rownames(GenoM)[SnpCallRate < 0.95] #Check for numbers of Samples below desired Call Rate
GenoM <- GenoM[, SnpCallRate > 0.60] # consider dropping samples below certain call rate if affecting data
hist(SNPCallRate, breaks = 20, col = "tan", main = "SNP Call Rate", xlab = "Call Rate")
hist(IndivCallRate, breaks = 20, col = "tan", main = "Individual Call Rate", xlab = "Call Rate")

#DUP/OH/ME; Test for different error values, ideally low values
CalcMaxMismatch(Err = 0.1, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #86/6/29
CalcMaxMismatch(Err = 0.01, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #4/0/0
CalcMaxMismatch(Err = 0.0001, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158)#0/0/0
CalcMaxMismatch(Err = 0.002, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #0/0/0
CalcMaxMismatch(Err = 0.003, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #0/0/0
CalcMaxMismatch(Err = 0.004, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #0/0/0
CalcMaxMismatch(Err = 0.005, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #1/0/0
CalcMaxMismatch(Err = 0.001, MAF = SnpStats(GenoFile, Plot = TRUE)[,"AF"], ErrFlavour = "version2.0", qntl = 0.005263158) #0/0/0

#Optional section, try different error values on simulated genotypes based on #SNP and #Individuals
#qntl = 0.005263158 from n = individuals checked by GenoConvert (95 for me) & Q = 0.5 with qntl = Q/N
Summary_Stats<- summary(SnpStats(GenoFile))
View(Summary_Stats) # My Data = AF Min: 0.3, AF Max: 0.5; Mis Min & Max: 0
MAFrange<-c(0.3, 0.5) # Based on AF Min & AF Max
Snp<-615
Indv<-120
Error_values1 <- c(0.01)  
Error_values2 <- c(0.001)  
Error_values3 <- c(0.0001)
Error_values4 <- c(0.00001)
Error_values5 <- c(0.000001) # You can test for desired error rate
Call_range <- c(0.8, 1) #Pick desired call rate (SNP present in % of individuals)
Missing_range <- c(0, 0) # Based on Mis Min & Mis Max

# Repeat for test 1 - 5, shown is test 5, using Error_values5
Simulated_genotype5 <- SimGeno(
  Pedigree = LizardPed,  # Replace with your actual pedigree data
  nSnp = Snp,
  MAF = runif(1, min = MAFrange[1], max = MAFrange[2]),
  SnpError = Error_values5,
  CallRate = runif(1, min = Call_range[1], max = Call_range[2]),
  ParMis = runif(1, min = Missing_range[1], max = Missing_range[2]),
  ReturnStats = TRUE  # Set ReturnStats to TRUE if you want statistics
)

# Repeat for test 1 - 5, shown is test 5, using Simulated_genotype5
Geno_matrix5 <- Simulated_genotype5$SGeno

# Repeat for test 1 - 5, shown is test 5, using Geno_matrix5 and Error_values5
Reconstructed_pedigree5 <- sequoia(
  GenoM = Geno_matrix5, Err = Error_values5,
  LifeHistData = LizardLH  # Provide life history data here
)
# End of optional section

#Run multiple pedigree building test using different error rates
Test001 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.001,
                   Plot = FALSE, quiet = "verbose", UseAge = "extra")
#parentage: assigned 57 dams and 45 sires to 95 individuals; full ped: assigned 66 dams and 66 sires to 95 + 10 individuals (real + dummy)
Test002 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.002,
                   Plot = FALSE, quiet = "verbose", UseAge = "extra")
#parentage: assigned 58 dams and 49 sires to 95 individuals; full ped: assigned 67 dams and 67 sires to 95 + 9 individuals (real + dummy)
Test003 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.003,
                   Plot = FALSE, quiet = "verbose", UseAge = "extra")
#parentage: assigned 58 dams and 50 sires to 95 individuals; full ped: assigned 72 dams and 73 sires to 95 + 13 individuals (real + dummy)
Test004 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.004,
                   Plot = FALSE, quiet = "verbose", UseAge = "extra")
#parentage: assigned 58 dams and 52 sires to 95 individuals; full ped: assigned 71 dams and 74 sires to 95 + 11 individuals (real + dummy)
Test005 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.005,
                   Plot = FALSE, quiet = "verbose", UseAge = "extra")
#parentage: assigned 58 dams and 52 sires to 95 individuals; full ped: assigned 71 dams and 74 sires to 95 + 11 individuals (real + dummy)

# Repeat for test 1 - 5, shown is test 5
sumry_Test005 <- SummarySeq(Test005, Plot=FALSE)
tmp5 <- apply(sumry_Test005$ParentCount['Genotyped',,,],
              MARGIN = c('parentSex', 'parentCat'), FUN = sum)
props5 <- sweep(tmp5, MARGIN='parentCat', STATS = rowSums(tmp5), FUN = '/')
1 - props5[,'Genotyped'] / (props5[,'Genotyped'] + props5[,'Dummy']) # Will produce results similar to example shown below

sumry_Test005
tmp5
props5
# Test 1: 0.1363636  dam, 0.2063492  sire
# Test 2: 0.1343284  dam, 0.1846154  sire
# Test 3: 0.1944444  dam, 0.2394366  sire
# Test 4: 0.1830986  dam, 0.2222222  sire
# Test 5: 0.1830986  dam, 0.2222222  sire

# Repeat for test 1 - 5, shown is test 5
Estconf_Test005 <- EstConf(Pedigree = Test005$Pedigree,
                           LifeHistData = LizardLH,
                           args.sim = list(nSnp = Snp,   # no. in actual data, or what-if
                                           SnpError = 0.005,  # best estimate, or what-if
                                           CallRate=0.99,     # from Call Rate determined from histogram of SNPCallRate
                                           ParMis=c(0.1830986, 0.2222222)),  # Copy from above results respective to each test
                           args.seq = list(Err=0.005, Module="ped"),  # as in real run
                           nSim = 20,   # try-out number of simulations, proper run >=20 (10 if huge pedigree)
                           nCores=1)

# Repeat for test 1 - 5, shown is test 5 to check for quality
Estconf_Test005$ConfProb
Estconf_Test005$PedErrors
1 - apply(Estconf_Test005$PedErrors, MARGIN=c('cat','parent'), FUN=sum, na.rm=TRUE)

#After testing for simulation and deciding on error rates which you are most satisfied for
TestTassignPed1 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.001, Tassign = 1.0,
                           Plot = FALSE, quiet = "verbose", UseAge = "extra")
#assigned 57 dams and 45 sires to 95 individuals; Full Pedigree - assigned 69 dams and 66 sires to 95 + 11 individuals (real + dummy)

TestTassignPed2 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.001, Tassign = 1.0,
                           Plot = FALSE, quiet = "verbose", UseAge = "yes")
#assigned 57 dams and 45 sires to 95 individuals; Full Pedigree - assigned 69 dams and 66 sires to 95 + 11 individuals (real + dummy)

TestTassignPed3 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.001, Tassign = 0.5,
                           Plot = FALSE, quiet = "verbose", UseAge = "extra")
#assigned 57 dams and 45 sires to 95 individuals; Full Pedigree - assigned 66 dams and 66 sires to 95 + 10 individuals (real + dummy)

TestTassignPed4 <- sequoia(GenoM = GenoFile, LifeHistData = LizardLH, Module = "ped", Err = 0.001, Tassign = 0.5,
                           Plot = FALSE, quiet = "verbose", UseAge = "yes")
#assigned 57 dams and 45 sires to 95 individuals; Full Pedigree - assigned 66 dams and 63 sires to 95 + 11 individuals (real + dummy)

#Compare Test pedigree to one another
TestTassignPedCompExtra <- PedCompare(TestTassignPed1$Pedigree, TestTassignPed3$Pedigree)
TestTassignPedCompExtra$Mismatch
TestTassignPedCompYes <- PedCompare(TestTassignPed2$Pedigree, TestTassignPed4$Pedigree)
TestTassignPedCompYes$Mismatch

#Compare test pedigree to known pedigree
Compare1.0extra<-PedCompare(TestTassignPed1$Pedigree, LizardPed)
Compare1.0yes<-PedCompare(TestTassignPed2$Pedigree, LizardPed)
Compare0.5extra<-PedCompare(TestTassignPed3$Pedigree, LizardPed)
Compare0.5yes<-PedCompare(TestTassignPed4$Pedigree, LizardPed)

