## This is psyllid 16S data analysis using dada2 package.
## all the conserved region primers were removed, low quality window trimmed by Trimmomatic
## the input are the reads that have gone through several rounds of trimming to size of 403 bp
## 403 bp was determined by the histogram of the merged reads done by Flash program

library(dada2)
library(vegan)
library(ggplot2)
library(magrittr)
library(reshape2)
library(phyloseq)

## Section I.  define file path and read data. 
wd <- "~/potato_psyllid_rad_tag/2016/16S_microbiome/dada2/"
setwd(wd)
#load("~/Back_up_win-2-25_2013/16S_psy/dada_files/taxa_rdp.rda")

ps <- readRDS("./dada_files/ps_silv.rds")


## reading new meta data
new_meta <- read.csv("./meta.csv", header = T)
adm <- read.csv("../../../Combined_data_12_16/Admixture/admixture_all_combined.txt",
                sep = "\t")
## read pca and admixture files

pca <- read.csv("../../../Combined_data_12_16/PCA/psy_2017_combined.evec",
                sep = "", header = T)

pca <- pca[, c(1:3)]
colnames(pca) <- c("sample_ID", "pc1", "pc2")
pca$sample_ID <- gsub("L5_", "L5-", pca$sample_ID)

new_meta$Pooling.ID <- gsub("_", "", new_meta$Pooling.ID)
mer <- merge(new_meta, adm, by.x= "INDV_SNP", by.y= "INDV", all.x =  T)
mer <- merge(mer, pca, by.x= "final_vcf_sample_names", by.y = "sample_ID", all.x=T)

rm(adm, pca)
mer <- mer[, c("Pooling.ID","final_vcf_sample_names",
               "k2_1", "k2_2", "k3_1", "k3_2", "k3_3", "pc1","pc2")]
## access the sample info
sam_info <- sample_data(ps@sam_data)

sam_info$k2_1 <- mer$k2_1[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$k2_2 <- mer$k2_2[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$k3_1 <- mer$k3_1[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$k3_2 <- mer$k3_2[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$k3_3 <- mer$k3_3[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$k2_major <- ifelse(sam_info$k2_1 >=0.5, "group1", "group2") %>%
                        as.factor()

sam_info$pca1 <- mer$pc1[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$pca2 <- mer$pc2[match(sam_info$Pooling.ID, mer$Pooling.ID)]
sam_info$SNP_ID <- mer$final_vcf_sample_names[match(sam_info$Pooling.ID,
                                                    mer$Pooling.ID)]

ps@sam_data <- sam_info


## Section II. next filter data based on abundance and contamination. 

ps_sub <- ps %>%
  subset_taxa(
    !is.na(Phylum)          & 
      Class != "Chloroplast" &
      !is.na(Class)   &
      Family != "Mitochondria"
  )
ps_sub

ps <- ps_sub
rm(ps_sub)

## change this cutoff value to e-3 or e-2. or whatever we want. 
filterfun1 = function(x){
  x[(x / sum(x)) < (1e-2)] <- 0
  return(x)
}
ps_transform <- transform_sample_counts(ps, fun = filterfun1)

ps_ra_fil <- filter_taxa(ps_transform , function(x) mean(x) > 0, TRUE)

### below is changed from "dada2_psy.R" 
filt <- genefilter_sample(ps_ra_fil, filterfun_sample(function(x) x >= 5), A=3)

ps2 <- prune_taxa(filt,ps_ra_fil)
ps2  
rm(ps_ra_fil, ps_transform, ps)

###### Section III, relative abundance and plotting
### outptu all the taxa and relative abundance ###### tutorial from 
### http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html ####
#### depends on what we want, we can adjust the ps here. we can use ps1 or ps 
#### if no relative abuandance is need. then use "ps2" and skip the below section

psr <- transform_sample_counts(ps2,function(x) x/sum(x))


psr3 <- subset_samples(psr, sample_names(psr)!= "Fu52")
psr3 <- subset_samples(psr3, sample_names(psr3)!= "Fu53")

taxa_names(psr3) <- paste0("taxon", seq(1:ntaxa(psr3)))

hap_order <- psr3@sam_data[order(psr3@sam_data$hap), ]$Pooling.ID
hap_order_snpid <-  psr3@sam_data[order(psr3@sam_data$hap), ]$SNP_ID
#host_order <- psr3@sam_data[order(psr3@sam_data$host), ]$Pooling.ID

## nightshade samples were 43,44,45,46,54,55,56,57,6
#psr@sam_data$num <- gsub("Fu","", psr@sam_data$Pooling.ID)
#psr@sam_data$num <- as.numeric(psr@sam_data$num)
#pooling_order <- psr@sam_data[order(psr@sam_data$num), ]$Pooling.ID

mdf = psmelt(psr3)
mdf = mdf[, c("OTU", "Sample", "Abundance")]
mdf$OTU <- as.factor(mdf$OTU)
mdf$Sample <- as.factor(mdf$Sample)

species_color <- c("lightblue2","mediumblue", "darkred", "maroon3",
                   "greenyellow", "seagreen2", "orange")
tax_order <- c("taxon1", "taxon7", "taxon5","taxon6", "taxon2", "taxon3",
               "taxon4")

p_sp <- ggplot(mdf, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color= "black") +
  scale_x_discrete(limits = hap_order, 
  #                 labels = c(rep("", 77))) +
                   labels = "") +
  scale_fill_manual(values = species_color, limits= tax_order, 
                    name = "") +
  labs(x = "", y = "Relative abundance") +
  theme_classic(base_size = 12) +
  theme(axis.text.x= element_text(angle = 90, hjust = 0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(rep(0.5, 4)), "cm"),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(angle=90, vjust=0.5, hjust = 1)) 
  
rm(tax, tax_order, p_sp, species_color, mdf, filt)



### this is the point to run PS object to all the filtering process to remove tax and sample 
## with different species of psyllids
## run the code above.... 

#### try to deal with the genetic grouping idea, updated on April 27,2018
psr3@sam_data$k2_major <- as.factor(psr3@sam_data$k2_major)

k2_order <- psr3@sam_data[order(psr3@sam_data$k2_major), ]$Pooling.ID
k2_order_snp <- sam_info$SNP_ID[match(k2_order, sam_info$Pooling.ID)]

mdf = psmelt(psr3)
mdf = mdf[, c("OTU", "Sample", "Abundance")]
mdf$OTU <- as.factor(mdf$OTU)
mdf$Sample <- as.factor(mdf$Sample)

species_color <- c("lightblue2","mediumblue", "darkred", "maroon3",
                   "greenyellow", "seagreen2", "orange")
tax_order <- c("taxon1", "taxon7", "taxon5","taxon6", "taxon2", "taxon3",
               "taxon4")

genetic <- ggplot(mdf, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", color= "black") +
  scale_x_discrete(limits = k2_order, 
               labels = c(rep("g1", 26), rep("g2", 46), rep("unknown","5"))) +
  #labels = host_order) +
  scale_fill_manual(values = species_color, limits= tax_order, 
                    name = "") +
  labs(x = "", y = "Relative abundance") +
  theme_classic(base_size = 12) +
  theme(axis.text.x= element_text(angle = 90, hjust = 0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(rep(0.5, 4)), "cm"),
        axis.ticks.x = element_blank(),
        axis.text.x =element_text(angle=90, vjust=0.5, hjust = 1)) 



## next plot the PCA based on the genetic grouping
gen_pca <- ggplot(sam_info, aes(x = pca1, y = pca2, color = k2_major)) +
  geom_point(size = 1.05, stroke = 1.02,alpha = 0.95) +
  scale_color_manual(values = c("firebrick4", "darkorange","gray75"),     
                     breaks = c("group1", "group2", ""),
                     labels=c("Genetic group 1", "Genetic group 2", 
                              ""))      +
  theme(legend.title=element_blank()) +
  labs(x = "PC 1 (15.37%)", y= "PC 2 (1.55%)") +
  #  theme(legend.text=element_text(size=12)) + 
  theme_classic(base_size= 12) +
  theme(legend.title=element_blank()) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) 
  
#### label the gen pca plot

gen_pca_lab <- ggplot(sam_info, aes(x = pca1, y = pca2, color = k2_major,
                                    label=SNP_ID)) +
  geom_point(size = 1.05, stroke = 1.02,alpha = 0.95) +
  geom_text(size = 2, check_overlap = T) +
  scale_color_manual(values = c("firebrick4", "darkorange","gray75"),     
                     breaks = c("group1", "group2", ""),
                     labels=c("Genetic group 1", "Genetic group 2", 
                              ""))      +
  theme(legend.title=element_blank()) +
  labs(x = "PC 1 (15.37%)", y= "PC 2 (1.55%)") +
  #  theme(legend.text=element_text(size=12)) + 
  theme_classic(base_size= 12) +
  theme(legend.title=element_blank()) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

        text = element_text(size = 4))
  theme(text = element_text(size = 10))


## combine with larger PCA plot
pca <- read.csv("../../../Combined_data_12_16/PCA/psy_2017_combined.evec",
                sep = "", header = T)
pca <- pca[, c(1:3)]
colnames(pca) <- c("sample_ID", "pc1", "pc2")
pca$sample_ID <- gsub("L5_", "L5-", pca$sample_ID)
pca_16 <- pca[c(1:508), ]
group1 <- sam_info[sam_info$k2_major =="group1", ]$SNP_ID
group2 <- sam_info[sam_info$k2_major =="group2", ]$SNP_ID
pca_16$color <- ifelse(pca_16$sample_ID %in% group1, "group1",
                       ifelse(pca_16$sample_ID %in% group2, "group2", "no")) 

pca_16$color <- as.factor(pca_16$color)
pca_16$size <- ifelse(pca_16$color=="group1" | pca_16$color=="group2",
                      "large", "small")

pca_all <- ggplot(pca_16, aes(x = pc1, y = pc2, color = color)) +
  geom_point(size = 1.02, stroke = 1.02, alpha = 0.85 ) + 
  scale_color_manual(values = c("firebrick4", "darkorange","gray75"),     
                     breaks = c("group1", "group2", "no"),
                    labels=c("Genetic group 1", "Genetic group 2", 
                              "NA for microbiome"))      +
  theme(legend.title=element_blank()) +
  labs(x = "PC 1 (15.37%)", y= "PC 2 (1.55%)") +
#  theme(legend.text=element_text(size=12)) + 
  theme_classic(base_size= 12) +
  theme(legend.title=element_blank()) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) 


## ###### Section IV. permutational amova 
## updated on Jul 5, 2018. PAMOVA with the new geneitc group
### next haplotype, first remove samples with unknown haplotypes ##
## can not use the relative abudnace phyloseq file. 
ps_gen <- subset_samples(ps2, sample_names(ps2)!= "Fu52")
ps_gen <- subset_samples(ps_gen, sample_names(ps_gen)!= "Fu53")
## there are a few 'NA' samples need to be removed
rem <- rownames(ps_gen@sam_data)[c(which(is.na(ps_gen@sam_data$k2_major)))]
ps_gen <- subset_samples(ps_gen, !sample_names(ps_gen) %in% rem)

ps_gen_scale <- scale_reads(ps_gen, min(sample_sums(ps_gen)))


ps_gen_bray <- phyloseq::distance(ps_gen_scale, method = "bray")
adonis(ps_gen_bray ~ k2_major, data = data.frame(sample_data(ps_gen_scale)))

beta_gen <- betadisper(ps_gen_bray, data.frame(sample_data(ps_gen))$k2_major)
permutest(beta_gen)
## the p-value is 0.003, so host is significant 


#### Section V: other functions used in the analysis ########
##### this is downloaded from https://github.com/michberr/MicrobeMiseq/blob/master/R/miseqR.R
# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}



