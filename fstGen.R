library(fst)
library(qtl2)
library(broman)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(qtlcharts)
library(data.table)

cat("Running Clean-o-type...\n")
cat("\nCalling script: fstGen.R\n")

#Run script from folder with the 3 input files 
#Form a list of filepaths where [[1]] == cross2, 
#[[2]] == codes, and [[3]] == Final Report
pathIn <- list(
  crossPath = list.files("../input/", pattern = "cross", full.names = T, ignore.case = T),
  codePath  = list.files("../input/", pattern = "code", full.names = T, ignore.case = T),
  reportPath= list.files("../input/", pattern = "report", full.names = T, ignore.case = T)
)


if(!all(sapply(pathIn, function(x) length(x) > 0))){
  cat("Error: One or more inputs could not be found.\n")
  print(sapply(pathIn, function(x) length(x) > 0))
}else{
  cat("\nFound the following inputs: \n")
  print(sapply(pathIn, function(x) length(x) > 0))
}

#Load cross
cross <- readRDS(pathIn[[1]])
#Load codes
codes <- read.csv(pathIn[[2]], skip = 3, stringsAsFactors = F) 
#Load a final report - read in the first 20 lines
tmp <- readLines(pathIn[[3]],n = 20)

skip = grep("^SNP Name\tSample",tmp) - 1
rm(tmp)

if(skip == 0){
  stop("Error: No headers detected within 20 lines.
\nCheck file Final Report and try again.")
}

cat("\nReading in Final Report. Please wait...\n")
finrep <- fread(pathIn[[3]], data.table = FALSE, skip = skip)
finrep <- finrep[finrep[[1]] %in% codes[[1]],]

cat("First 5 Sample IDs...\n")
viewID <- head(unique(finrep[[2]]), n = 5)
cat(viewID, sep = "\n")

altIDs <- "empty"
while(!altIDs %in% c("Y","N")){
  altIDs <- toupper(
    readline("Would you like to alter sample IDs(Y/N)? "))
}

if(altIDs == "Y"){
  desig <- readline("Prefix sample IDs with...: ")
  samp <- paste0(desig, viewID)
  cat("IDs will now look like: \n")
  cat(samp, sep = "\n")
  conf1 <- toupper(readline("Confirm change (Y/N): "))
  
  if(conf1 == "Y"){
    finrep[[2]] <- paste0(desig,finrep[[2]])
    cat("IDs have been changed.\n")
  }else{
    cat("Aborting...\n
        Please try again.")
  }
}

cat("Pulling Intensities...\n")

# create matrices that are snps x samples
snps <- unique(finrep[,"SNP Name"])
samples <- unique(finrep[,"Sample ID"])
X <- Y <- matrix(ncol=length(samples), nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
  message(i, " of ", length(samples))
  tmp <- finrep[finrep[,"Sample ID"]==samples[i],]
  X[,samples[i]] <- tmp[,"X"]
  Y[,samples[i]] <- tmp[,"Y"]
}
cat("\nFinished. \n")
# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps), seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

# write to fst file, maximally compressed
cat("Writing out FST...")
cat("Finished.\n")
write.fst(result, "../output/objects/DOWL_intensities.fst", compress=100)
rm(finrep)

#SNPs on X chromosome
Xsnp <- codes$marker[codes$chr == "X"]
Ysnp <- codes$marker[codes$chr == "Y"]

cat("\nGenerating Sex Intensity Diagnostic...\n")
#Process X Intensity
Xres <- result[result$snp %in% Xsnp,]
#T-test markers for significant intensity by sex
Xint <- (Xres[Xres$channel == "X", -c(1,2)] + Xres[Xres$channel == "Y",-c(1,2)])/2
rownames(Xint) <- Xres$snp[Xres$channel == "X"]

xpval <- apply(t(Xint), 2, function(a) t.test(a~cross$covar[match(colnames(Xint),
                                                                  rownames(cross$covar)),
                                                            "sex"], 
                                              na.rm = T)$p.value)
Xavg <- apply(Xint[xpval <= 0.05,], 2, mean)

#Process Y Intensity
Yres <- result[result$snp %in% Ysnp,]
rm(result)
Yint <- (Yres[Yres$channel == "X", -c(1,2)] + Yres[Yres$channel == "Y",-c(1,2)])/2
rownames(Yint) <- Yres$snp[Yres$channel == "Y"]
#T-test for significant intensity by sex
ypval <- apply(t(Yint), 2, function(a) t.test(a~cross$covar[match(colnames(Yint),
                                                                  rownames(cross$covar)),
                                                            "sex"], 
                                              na.rm = T)$p.value)
Yavg <- apply(Yint[ypval <= 0.05,], 2, mean)

intDF <- data.frame(Xavg = Xavg, Yavg = Yavg)
#Add in sex 
intDF$sex <- cross$covar[match(rownames(intDF),rownames(cross$covar)),"sex"]
intDF <- intDF[!is.na(intDF$sex),]
phetX <- rowSums(cross$geno$X == 2)/rowSums(cross$geno$X != 0)
intDF$phetX <- phetX[match(rownames(intDF),names(phetX))]

flag <- character()
for(i in unique(intDF$sex)){
  tmp <- intDF[intDF$sex == i,]
  xthresh_hi <- mean(tmp$Xavg) + (3*sd(tmp$Xavg))
  xthresh_lo <- mean(tmp$Xavg) - (3*sd(tmp$Xavg))
  ythresh_hi <- mean(tmp$Yavg) + (3*sd(tmp$Yavg))
  ythresh_lo <- mean(tmp$Yavg) - (3*sd(tmp$Yavg))

  flag <- unique(c(flag, rownames(tmp[which(tmp$Xavg >= xthresh_hi |
              tmp$Xavg <= xthresh_lo |
              tmp$Yavg >= ythresh_hi |
              tmp$Yavg <= ythresh_lo),])))

}


a <- ggplot(intDF, aes(x = Xavg, y = Yavg, color = sex))+ 
  geom_point()+
  geom_label_repel(data = intDF[rownames(intDF) %in% flag,],
                   aes(label = rownames(intDF[rownames(intDF) %in% flag,])),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   show.legend = F)+
  labs(title = "Sex Chromosome Diagnostic",
       x = "Average X Intensity", 
       y = "Average Y Intensity")+
  theme_bw()

#When you get back - add a flagger for ggrepel here - just swap out y for prop Het X
flag2 <- character()
for(i in unique(intDF$sex)){
  tmp <- intDF[intDF$sex == i,]
  xthresh_hi <- mean(tmp$Xavg) + (3*sd(tmp$Xavg))
  xthresh_lo <- mean(tmp$Xavg) - (3*sd(tmp$Xavg))
  phetXthresh_hi <- mean(tmp$phetX) + (3*sd(tmp$phetX))
  phetXthresh_lo <- mean(tmp$phetX) - (3*sd(tmp$phetX))
  
  flag2 <- unique(c(flag2, rownames(tmp[which((tmp$Xavg >= xthresh_hi |
                                        tmp$Xavg <= xthresh_lo) | (
                                        tmp$phetX >= phetXthresh_hi |
                                        tmp$phetX <= phetXthresh_lo)),])))
  
}



b <- ggplot(intDF, aes(x = Xavg, y = phetX, color = sex))+ 
  geom_point()+
  labs(x = "Average X Intensity", 
       y = "Proportion of Heterogeneity on X Chr")+
  geom_label_repel(data = intDF[rownames(intDF) %in% flag2,],
                   aes(label = rownames(intDF[rownames(intDF) %in% flag2,])),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   show.legend = F)+
  theme_bw()




ggarrange(a,b, 
          ncol = 1, nrow = 2, 
          common.legend = T,
          legend = "right") %>% ggexport(filename = "../output/plots/SexDiag.pdf")


invisible(gc())
cat("fstGen.R: Status: Complete")













