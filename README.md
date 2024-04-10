Codes used to create filtered SNPs, Pedigree, Kinship Matrix, and Sequoia pedigree
Codes used for Sequoia Pedigree is adopted from Cassandra Rodriguez of the UC Davis - Mammalian Ecology Conservation Unit\

My Steps:
1. Create pedigree figure using kinship2 based on information provided by the zoo. Inbreeding is determined from relationship and pairings.
2. Filter SNPs using plink, make kinship matrix of founder using plink2-KING
3. Make genomic pedigree using sequoia
4. Check kinship matrix and inferred pedigree of founders for matches
