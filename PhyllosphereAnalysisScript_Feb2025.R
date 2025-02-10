
#### Install all appropriate packages ####
install.packages("rlang")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("remotes")
install.packages("writexl")
install.packages("ape")
install.packages("vegan")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("paletteer")

remotes::install_github("vmikk/metagMisc")

### Takes a while to install the following package 
remotes::install_github("HuaZou/MicrobiomeAnalysis", force = T)

BiocManager::install("dada2", force = TRUE)
BiocManager::install("phyloseq", force = TRUE)
BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("phangorn", force = TRUE)
BiocManager::install("DECIPHER", force = TRUE)

require(writexl);require(dada2);require(phyloseq); require(Biostrings); require(vegan);
require(ape);require(phangorn);require(DECIPHER);require(vegan);
require(tidyverse); require(MicrobiomeAnalysis); require(metagMisc); require(RColorBrewer);
require(paletteer)

## Run this code before creating graphs!
theme_set(theme_bw())

#### Microbiome Data is the folder with all of mt fastq files for each sample.####
path<-"/Users/cdice/Microbiome Data/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2.fastq"))
names(filtFs) <- sample.names
names(filtFs)
names(filtRs) <- sample.names


# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230),
                     
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     
                     compress=TRUE, multithread=FALSE, trimLeft=21)

head(out)
out

errF <- learnErrors(filtFs, multithread=TRUE)
errF
errR <- learnErrors(filtRs, multithread=TRUE)
errR
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
mergers

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
view(track)

trackdata<-as.data.frame(track)
mean(trackdata$nonchim)
trackdata$inoculantstatus<-as.factor(c("N","N","N","N","N","N","N","Y","Y","Y","Y","Y","Y","Y","Y"))

summary(aov(trackdata$nonchim~trackdata$inoculantstatus))

#### Inoculated leaves have a significantly higher number of reads, so ratification is recommended ####

#####BEGINNING OF PHYLOGENETIC TREE CODE- code used to create the phylogeny tree of taxa. This is NOT NEEDED to run the dada 2####
##### pipeline and can be skipped. I had 21 samples, and the final command takes about 8 hours on my laptop so be warned####

###Not sure what the following two lines of code were used for - I would ignore if you want to use the rest of the code
write.csv(track, "track.csv")
sum(track$nonchim)

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

fitGTR

#####END OF PHYLOGENETIC TREE CODE#######

###I previously had some issues calling in the Silva Database but this should be the correct format####
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa")


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print


### This code was used to denote my two treatments

samples.out <- rownames(seqtab.nochim)
samples.out
subject <- sapply(samples.out, `[`, 1)
subject
Inoculation<-substr(subject,3,3)
Inoculation
subject <- substr(subject,3,4)
subject
samdf <- data.frame(Subject=subject, Inoculation=Inoculation) %>% 
  mutate(Color = case_when(Inoculation == "C" ~ "red", Inoculation == "D" ~ "red", Inoculation == "A" ~ "blue", Inoculation == "B" ~ "blue",), 
         Inoculation = case_when(Inoculation == "C" ~ "Inoculated", Inoculation == "D"~"Inoculated", TRUE ~ "Control"))
rownames(samdf) <- samples.out
rownames(samdf)

samdf

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa), phy_tree(fitGTR$tree))

tax_table(ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps <- subset_taxa(ps, Order != "Chloroplast" | is.na(Order))
ps <- subset_taxa(ps, Family != "Mitochondria" | is.na(Family))

mean(sample_sums(ps))

ps

readcount<-as.data.frame(sample_sums(ps))
mean(readcount[1:7,])
mean(readcount[8:15,])

t.test(readcount[1:7,], (readcount[8:15,]))

### Rarefaction of datasets with or without methylobacterium 
ps_meth<-ps
ps_meth

ps_meth<-rarefy_even_depth(ps_meth, sample.size = min(sample_sums(ps_meth)),
                      rngseed = 101010, replace = F, trimOTUs = F, verbose = TRUE)

sample_sums(ps_meth)

ps_no_meth<-subset_taxa(ps, Genus != "Methylobacterium-Methylorubrum" | is.na(Genus))

sample_sums(ps_no_meth)
 
readcount_no_meth<-as.data.frame(sample_sums(ps_no_meth))
mean(readcount_no_meth[1:7,])
mean(readcount_no_meth[8:15,])
mean(readcount_no_meth[,])

rotate_meth<-t(otu_table(ps_meth))

ASVTable_meth<-cbind(rotate_meth, tax_table(ps_meth))
df_meth<-ASVTable_meth
df
write.csv(df_meth, "DADA2_ASVTableRound2_meth.csv")

rotate_no_meth<-t(otu_table(ps_no_meth))

ASVTable_no_meth<-cbind(rotate_no_meth, tax_table(ps_no_meth))
df_no_meth<-ASVTable_no_meth
df
write.csv(df_no_meth, "DADA2_ASVTableRound2_no_meth.csv")

library(tidyverse)
library(vegan)

wUF.ordu_no_meth = ordinate(ps_no_meth, method="NMDS", distance="unifrac", weighted=TRUE)
stressplot(wUF.ordu_no_meth)
wUF.ordu_no_meth$stress

bray.ordu_no_meth = ordinate(ps_no_meth, method="NMDS", distance="bray", weighted=FALSE)
jaccard.ordu_no_meth = ordinate(ps_no_meth, method="NMDS", distance="jaccard", weighted=FALSE)

plot_ordination(ps_no_meth, wUF.ordu_no_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  

plot_ordination(ps_no_meth, bray.ordu_no_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  

plot_ordination(ps_no_meth, jaccard.ordu_no_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  


###Ordinate - with methylobacterium 

wUF.ordu_meth = ordinate(ps_meth, method="NMDS", distance="unifrac", weighted=TRUE)
stressplot(wUF.ordu_meth)
wUF.ordu_meth$stress

bray.ordu_meth = ordinate(ps_meth, method="NMDS", distance="bray", weighted=FALSE)
jaccard.ordu_meth = ordinate(ps_meth, method="NMDS", distance="jaccard", weighted=FALSE)

plot_ordination(ps_meth, wUF.ordu_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  

plot_ordination(ps_meth, bray.ordu_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  

plot_ordination(ps_meth, jaccard.ordu_meth, type="sites", color="Inoculation")  + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  stat_ellipse() + scale_color_manual(values = c("red", "blue"), labels=c('Control', 'Inoculated'))  


####alpha diversity tests

diversity_no_meth<-estimate_richness(ps_no_meth, split = TRUE)
write.csv(diversity_no_meth, "diversitymetricsRound2_no_meth.csv")

diversity_meth<-estimate_richness(ps_meth, split = TRUE)
write.csv(diversity_meth, "diversitymetricsRound2_meth.csv")

div_no_meth<-read.csv("diversitymetricsRound2_no_meth.csv", header = TRUE)
div_meth<-read.csv("diversitymetricsRound2_meth.csv", header = TRUE)

summary(Shanmod1<-aov(Shannon~Inoculation, data = div_no_meth))
summary(Simpmod1<-aov(Simpson~Inoculation, data = div_no_meth))

summary(Shanmod1<-aov(Shannon~Inoculation, data = div_meth))
summary(Simpmod1<-aov(Simpson~Inoculation, data = div_meth))

### Beta Diversity tests

library(MicrobiomeAnalysis)
packageVersion("vegan")

anosim(otu_table(ps_no_meth),group = Inoculation, distance = "bray")
anosim(otu_table(ps_no_meth),group = Inoculation, distance = "jaccard")

anosim(otu_table(ps_meth),group = Inoculation, distance = "bray")
anosim(otu_table(ps_meth),group = Inoculation, distance = "jaccard")

run_ANOSIM(ps_no_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")
run_ANOSIM(ps_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")

run_ANOSIM(ps_no_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")
run_ANOSIM(ps_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")

run_ANOSIM(ps_no_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")
run_ANOSIM(ps_meth, level = "Genus", variable = "Inoculation", method = "wunifrac")

#### Relative Abundance Graph - No Methylobacteria 

ggtheme <- theme_bw() + theme(panel.grid = element_blank())

leafotus <- read.csv("DADA2_ASVTableRound2_no_meth.csv")
view(leafotus)

leafotus_alldat_tmp <- leafotus %>%
  select(Genus, R2A1:R2D6) %>%
  group_by(Genus) %>%
  summarize_at(vars(R2A1:R2D6), sum) %>%
  mutate(Genus = case_when(is.na(Genus) == TRUE ~ "Other", TRUE ~ Genus)) %>%
  rowwise() %>% mutate(rowsum = sum(c_across(R2A1:R2D6))) %>%
  arrange(desc(rowsum)) %>%
  select(!rowsum)

view(leafotus_alldat_tmp)

leafotus_other <- as.data.frame(c(Genus = "Other",
                            colSums(leafotus_alldat_tmp[16:nrow(leafotus_alldat_tmp), 2:16]))) %>% t

leafotus_spec <- tibble(rbind(leafotus_alldat_tmp[1:15,], leafotus_other)) %>%
  mutate_at(vars(R2A1:R2D6), as.numeric) %>%
  pivot_longer(cols = R2A1:R2D6, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Treatment = case_when(str_detect(Sample, "[A|B]") ~ "Control", TRUE ~ "Inoculated")) %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Other", after = Inf))

leafotus_no_meth_rel<-aggregate(leafotus_spec$Abundance, list(leafotus_spec$Genus, leafotus_spec$Treatment), sum)

write.csv(leafotus_no_meth_rel, "RelativeAbundanceDataRound2_no_meth_top15.csv")

ggplot(data = leafotus_spec, aes(fill=Genus, y=Abundance, x=Treatment)) +
  ggtheme + labs(y= "Relative Abundance") +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#FFCC00FF","#FF0066FF", "#328C97FF", "#D1AAC2FF", "#A5506DFF", "#B3E0BFFF", "#2A9D3DFF", 
                               "#EDF181FF", "#DB7003FF", "#FBA600FF", "#F8C1A6FF", "#A30000FF", "#FF3200FF",
                               "#011A51FF", "#8090B0FF"))
##### With Methylobacteria 

leafotus_meth <- read.csv("DADA2_ASVTableRound2_meth.csv")
view(leafotus_meth)

leafotus_alldat_tmp_meth <- leafotus_meth %>%
  select(Genus, R2A1:R2D6) %>%
  group_by(Genus) %>%
  summarize_at(vars(R2A1:R2D6), sum) %>%
  mutate(Genus = case_when(is.na(Genus) == TRUE ~ "Other", TRUE ~ Genus)) %>%
  rowwise() %>% mutate(rowsum = sum(c_across(R2A1:R2D6))) %>%
  arrange(desc(rowsum)) %>%
  select(!rowsum)

view(leafotus_alldat_tmp_meth)

leafotus_other_meth <- as.data.frame(c(Genus = "Other",
                                  colSums(leafotus_alldat_tmp_meth[16:nrow(leafotus_alldat_tmp_meth), 2:16]))) %>% t


leafotus_spec_meth <- tibble(rbind(leafotus_alldat_tmp_meth[1:15,], leafotus_other)) %>%
  mutate_at(vars(R2A1:R2D6), as.numeric) %>%
  pivot_longer(cols = R2A1:R2D6, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Treatment = case_when(str_detect(Sample, "[A|B]") ~ "Control", TRUE ~ "Inoculated")) %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Other", after = Inf))

leafotus_meth_rel<-aggregate(leafotus_spec_meth$Abundance, list(leafotus_spec_meth$Genus, leafotus_spec_meth$Treatment), sum)

write.csv(leafotus_meth_rel, "RelativeAbundanceDataRound2_meth_top15.csv")

ggplot(data = leafotus_spec_meth, aes(fill=Genus, y=Abundance, x=Treatment)) +
  ggtheme + labs(y= "Relative Abundance") +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c(paletteer::paletteer_d("MoMAColors::Warhol"),
                               paletteer::paletteer_d("palettetown::ralts")))


ggplot(data = leafotus_spec_meth, aes(fill=Genus, y=Abundance, x=Treatment)) +
  ggtheme + labs(y= "Relative Abundance") +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#FF0066FF", "#328C97FF", "#D1AAC2FF", "blue","#A5506DFF", "#B3E0BFFF", "#2A9D3DFF", 
                               "#EDF181FF", "#EF5FAFFF", "#DB7003FF", "#FBA600FF", "#F8C1A6FF", "#FF3200FF",
                               "#011A51FF", "#8090B0FF"))


