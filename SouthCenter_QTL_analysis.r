#QTL analysis for SC F2 panel
#Matthew Wersebe, University of Oklahoma
#June 10, 2021
############################################################################
library(qtl2)
setwd("/home/weider/SC_QTL_mapping")

#read in the cross information from the control file:

SC_cross <- read_cross2("SouthCenter.yaml")

#Insert the psudeomarkers:
map <- insert_pseudomarkers(SC_cross$gmap, step = 1)

#Calculate genotype probabilities:

pr <- calc_genoprob(SC_cross, map, error_prob = 0.002, cores = 4)

#Calacule the allele probalities:

apr <- genoprob_to_alleleprob(pr)

#Genome Scan:

out1 <- scan1(pr, SC_cross$pheno, cores = 4)

#sanity check:
head(out1)

#Plot the LOD scores and peaks:
ymax <- maxlod(out1)

#PUE
plot(out1, map, lodcolumn=1, col = "slateblue", ylim = c(0, ymax*1.02))
legend("topleft", lwd = 2, col = "slateblue", colnames(out1), bg = "gray90")

#GR_DIff
plot(out1, map, lodcolumn=2, col = "darkblue", ylim = c(0, ymax*1.02))
legend("topright", lwd = 2, col = "darkblue", colnames(out1), bg = "gray90")

#Body_P
plot(out1, map, lodcolumn=3, col = "darkred", ylim = c(0, ymax*1.02))
legend("topleft", lwd = 2, col = "darkred", colnames(out1), bg = "gray90")



#look at the peaks:
peaks <- find_peaks(out1, map, threshold = 1.5, drop = 0.5)
dim(peaks)
peaks

#perform permutation test:

out1perm <- scan1perm(pr, SC_cross$pheno, cores = 4, n_perm = 1000)

hist(out1perm[,'PUE'], breaks = 50, xlab = "LOD", main = "LOD scores for PUE scan with threshold in red")
abline(v = summary(out1perm)[,'PUE'], col = 'red', lwd = 2)

hist(out1perm[,'GR_Diff'], breaks = 50, xlab = "LOD", main = "LOD scores for Growht Rate scan with threshold in red")
abline(v = summary(out1perm)[,'GR_Diff'], col = 'red', lwd = 2)


hist(out1perm[,'Body_P'], breaks = 50, xlab = "LOD", main = "LOD scores for Body %P scan with threshold in red")
abline(v = summary(out1perm)[,'Body_P'], col = 'red', lwd = 2)
summary(out1perm)
