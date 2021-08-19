###EdgeR
library(edgeR)
## Load Files individually into a data.frame
## Alternatively, you can merge the gene counts files manually or using shell scripting if you prefer that method and load them in to R directly
setwd("/Users/alexwaugh/Desktop/Hunt_Lab/JH_Repro_Ontogeny/ReadsPerGene_out")
temp = list.files(pattern="*.ReadsPerGene.out.tab")
temp2 <- sapply(strsplit(temp, split='.', fixed=TRUE), function(x) (x[1]))
myfiles = lapply(temp, read.delim,header=F)
myfiles_new <- myfiles

## merge first two columns into identifier column
for (i in 1:20){
  colnames(myfiles_new[[i]])[2] <- paste(temp2[i],"_counts",sep="")
  myfiles_new[[i]] <- myfiles_new[[i]][,c(1,2)]
}

## merge all data.frames by ID and reformat into gene_data.frame
merged <- Reduce(function(...) merge(..., by='V1', all.x=FALSE), myfiles_new)
row.names(merged) <- merged$V1

colnames(merged) <- c("V1",
                      #SRR7209532
                      "Sinv_F3",
                      #SRR7209533
                      "Sinv_G3",
                      #SRR7209534
                      "Sinv_F1",
                      #SRR7209535
                      "Sinv_N4",
                      #SRR7209536
                      "Sinv_F4",
                      #SRR7209537
                      "Sinv_N5",
                      #SRR7209538
                      "Sinv_F5",
                      #SRR7209539
                      "Sinv_N3",
                      #SRR7209540
                      "Sinv_Q3",
                      #SRR7209541
                      "Sinv_Q2",
                      #SRR7209542
                      "Sinv_Q1",
                      #SRR7209543
                      "Sinv_G5",
                      #SRR7209544
                      "Sinv_N2",
                      #SRR7209545
                      "Sinv_N1",
                      #SRR7209546
                      "Sinv_G2",
                      #SRR7209547
                      "Sinv_G1",
                      #SRR7209548
                      "Sinv_G4",
                      #SRR7209549
                      "Sinv_F2",
                      #SRR7209550
                      "Sinv_Q5",
                      #SRR7209551
                      "Sinv_Q4")
colnames(merged)
merged_ordered <- merged[,c(15,14,9,5,7,4,19,2,6,8,17,16,3,18,13,12,11,10,21,20)]
colnames(merged_ordered)
# create advanced group file that gives EdgeR info about your exp design
adv_group_Brain <- read.delim("/Users/alexwaugh/Desktop/Hunt_Lab/JH_Repro_Ontogeny/JH_Brain_adv_group.csv", row.names = 1, sep = ",")
# caste_ID
# specific_ID
# queen_class_ID

## Pairwise Comparisons
# subset brain and ovary data separately
y_caste <- DGEList(merged_ordered, group = adv_group_Brain$caste)


# build the pairwise design for caste
design_pair_caste <- model.matrix(~0+adv_group_Brain$caste, data=y_caste$Sample)

# filter genes with cpm < 1 in 6 of the 16 libraries in each comparison
keep_caste <- rowSums(cpm(y_caste) > 1) >= 5
y_caste <- y_caste[keep_caste, , keep.lib.sizes=FALSE]
table(keep_caste)


## Normalize samples by library size estimate various dispersions for EdgeR's statistical tests
#Caste
y_caste <- calcNormFactors(y_caste)
Background.genes_caste <- row.names(y_caste$counts)
y_caste <- estimateDisp(y_caste,design_pair_caste)
y_caste <- estimateGLMCommonDisp(y_caste, design_pair_caste)
y_caste <- estimateGLMTrendedDisp(y_caste, design_pair_caste)
y_caste <- estimateGLMTagwiseDisp(y_caste, design_pair_caste)
fit_caste <- glmQLFit(y_caste,design_pair_caste)


##Make Contrasts
colnames(design_pair_caste) <- c("queen", "worker")
my.contrasts_caste <- makeContrasts(
  queen_vs_worker=queen-worker,
  levels=design_pair_caste)


##### GLM QLF TEST ######
# Brain
glmQLFTest(fit_caste, contrast=my.contrasts_caste[,"queen_vs_worker"]) -> queen_vs_worker.test

########### topTags ###########
# Brain
topTags(glmQLFTest(fit_caste, contrast=my.contrasts_caste[,"queen_vs_worker"]), p.value = 1, n = 100000)$table -> QLF_queen_vs_worker.test

#write a usable table with cpm
y_caste_cpm <- cpm(y_caste)


######################################################

### DE Analysis

## isolate FDR corrected pvalues < 0.10, 0.05, 0.01, 0.001

###########Brain###########
subset(QLF_queen_vs_worker.test, FDR < 0.10) -> queen_vs_worker_10
subset(QLF_queen_vs_worker.test, FDR < 0.05) -> queen_vs_worker_05
subset(QLF_queen_vs_worker.test, FDR < 0.01) -> queen_vs_worker_01
subset(QLF_queen_vs_worker.test, FDR < 0.001) -> queen_vs_worker_001

########
