


# *******Script for data analyses ******************************** -----
# Manuscript:   Microbial Communities Associated with Cultivated Morel (Morchella spp.) Mushrooms
# Authors:      Longley R, Benucci GMN, ... Bonito G.
# Affiliation:  Michigan State University
# Journal:      FEMS Microbiology Ecology 
# Date:         January 23, 2019
# ******************************************************************** -----




#The original miseq dataset for this script contains data for two different datasets. The first few steps 
#are to remove contaminants and check for tag switching and are done on the entire dataset. Following these steps
#all other steps and figures are generated for the manuscript listed above only. 



# WORKING ENVIRONMENT SETUP ------------------------------------------------------------
options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
options(verbose=TRUE)

# >>> importing datasets --------------------------------------------------------------
library(ggplot2)
library(phyloseq)
library(vegan)
library(dplyr)
library(Biostrings)
library(ape)


# Import data ------------
#Fungi
otu_table_fungi <- read.csv("Suppl_file1_otu_table_fungi.csv", header=T, row.names =1)
otu_table_fungi_phy <-otu_table(otu_table_fungi,taxa_are_rows = TRUE)

tax_table_fungi <- read.csv("Suppl_file2_tax_table_fungi.csv", header=T, row.names =1)
tax_table_fungi_phy <- tax_table(as.matrix(tax_table_fungi))

sample_data_fungi <- read.csv("Suppl_file3_sample_data_fungi.csv", header=T, row.names =1)
sample_data_fungi_phy <-sample_data(sample_data_fungi)

otus_fungi <- readDNAStringSet("Suppl_file4_otus_fungi.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# create a phyloseq object for fungi
physeq_fungi <- merge_phyloseq(otu_table_fungi_phy, 
                               tax_table_fungi_phy,
                               sample_data_fungi_phy,
                               otus_fungi)

physeq_fungi
sample_data(physeq_fungi)
otu_table(physeq_fungi)
tax_table(physeq_fungi)

#import prokaryote data
otu_table_prokaryote <- read.csv("Suppl_file5_otu_table_prokaryote.csv", header=T, row.names =1)
otu_table_prokaryote_phy <-otu_table(otu_table_prokaryote,taxa_are_rows = TRUE)

tax_table_prokaryote <- read.csv("Suppl_file6_tax_table_prokaryote.csv", header=T, row.names =1)
tax_table_prokaryote_phy <- tax_table(as.matrix(tax_table_prokaryote))

sample_data_prokaryote <- read.csv("Suppl_file7_sample_data_prokaryote.csv", header=T, row.names =1)
sample_data_prokaryote_phy <-sample_data(sample_data_prokaryote)

otus_prokaryote <- readDNAStringSet("Suppl_file8_otus_prokaryote.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# create a phyloseq object for prokaryotes
physeq_prokaryote <- merge_phyloseq(otu_table_prokaryote_phy, 
                                    tax_table_prokaryote_phy,
                                    sample_data_prokaryote_phy,
                                    otus_prokaryote)

physeq_prokaryote
sample_data(physeq_prokaryote)
otu_table(physeq_prokaryote)
tax_table(physeq_prokaryote)


# extracting MOCK positive samples ---------------
physeq_fungi -> biom_ITS_uparse

biom_ITS_mock <- subset_samples(biom_ITS_uparse, Description%in%c("MOCK4"))
otu_table(biom_ITS_mock) <- otu_table(biom_ITS_mock)[which(rowSums(otu_table(biom_ITS_mock)) >= 1),] 
biom_ITS_mock
otu_table(biom_ITS_mock)
tax_table(biom_ITS_mock)
write.table(refseq(biom_ITS_mock), "MOCK1_Morel_ITS.fasta")

physeq_prokaryote -> biom_16s_uparse

biom_16s_mock <- subset_samples(biom_16s_uparse, Description%in%c("MOCK1"))
otu_table(biom_16s_mock) <- otu_table(biom_16s_mock)[which(rowSums(otu_table(biom_16s_mock)) >= 1),] 
biom_16s_mock
otu_table(biom_16s_mock)
tax_table(biom_16s_mock)
write.table(refseq(biom_16s_mock), "MOCK1_Morel_16s.fasta")


# >>> Filtering out OTUs ------------------------------
# Lindhal et al. 2013, tag switching - that's a good  one!'''
# Barberan et al. 2012, removing OTUs that appear in less than x samples


biom_ITS_uparse -> biom_ITS_uparse_filt
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 10),] ### PCR Errors 
biom_ITS_uparse_filt

biom_16s_uparse -> biom_16s_uparse_filt
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 10),] ### PCR Errors 
biom_16s_uparse_filt


# function to remove bad taxa -----------
# this will remove taxa from the dataset that appeared to be contaminants(as described in the manuscript)
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# remove ITS contaminants
badTaxa_ITS =c("OTU_17", "OTU_51")

biom_ITS_uparse_filt = remove_taxa(biom_ITS_uparse_filt, badTaxa_ITS)
biom_ITS_uparse_filt
sample_data(biom_ITS_uparse_filt)

# remove 16S contaminants
badTaxa_16s =c("OTU_260","OTU_1933", "OTU_245","OTU_1155","OTU_2789","OTU_2693","OTU_9938","OTU_2645",
               "OTU_348","OTU_1790","OTU_414","OTU_2021","OTU_11","OTU_3102","OTU_5590","OTU_1530",
               "OTU_111","OTU_811","OTU_78","OTU_9041","OTU_104","OTU_8262","OTU_2701","OTU_5","OTU_150",
               "OTU_112","OTU_68","OTU_460","OTU_161","OTU_786","OTU_883","OTU_532","OTU_1149","OTU_889",
               "OTU_1244","OTU_482","OTU_1001","OTU_1172","OTU_692","OTU_957","OTU_1127","OTU_4226",
               "OTU_601","OTU_977")

biom_16s_uparse_filt = remove_taxa(biom_16s_uparse_filt, badTaxa_16s)
biom_16s_uparse_filt
sample_data(biom_16s_uparse_filt)

# splitting out the indoor dataset used in the manuscript-----------------------------------

# Morel ITS indoor 
biom_ITS_uparse_in <- subset_samples(biom_ITS_uparse_filt, Study%in%c("Indoor"))
otu_table(biom_ITS_uparse_in) <- otu_table(biom_ITS_uparse_in)[which(rowSums(otu_table(biom_ITS_uparse_in)) >= 1),] 
biom_ITS_uparse_in
sample_data(biom_ITS_uparse_in)
otu_table(biom_ITS_uparse_in)

#Morel 16S Indoor
biom_16s_uparse_in <- subset_samples(biom_16s_uparse_filt, Study%in%c("Indoor"))
otu_table(biom_16s_uparse_in) <- otu_table(biom_16s_uparse_in)[which(rowSums(otu_table(biom_16s_uparse_in)) >= 1),] 
biom_16s_uparse_in
sample_data(biom_16s_uparse_in)
otu_table(biom_16s_uparse_in)


#stacked barplots----------

#prokaryotes

biom_16s_uparse_in_class <- biom_16s_uparse_in %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)
library(dplyr)

# Plot prokaryote bars
bar_prok = ggplot(biom_16s_uparse_in_class, aes(x = Label, y = Abundance, fill = Class)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits=c("Primordia 1", "Primordia 2", "Primordia 3", "Differentiation 1", "Differentiation 2","Differentiation 3", "Fundaments 1","Fundaments 2", "Fundaments 3", "Late fruiting 1", "Late fruiting 2", "Late fruiting 3", "Late non fruiting 1", "Late non fruiting 2", "Late non fruiting 3"))+
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Classes > 2%) \n") 
plot(bar_prok)

#fungi

biom_ITS_uparse_in_class <- biom_ITS_uparse_in %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by class


phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)



# Plot fungi bars
bar_fungi = ggplot(biom_ITS_uparse_in_class, aes(x = Label, y = Abundance, fill = Class)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits=c("Primordia 1", "Primordia 2", "Primordia 3", "Differentiation 1", "Differentiation 2","Differentiation 3", "Fundaments 1","Fundaments 2", "Fundaments 3", "Late fruiting 1", "Late fruiting 2", "Late fruiting 3", "Late non fruiting 1", "Late non fruiting 2", "Late non fruiting 3"))+
  scale_fill_manual(values = phylum_colors) +
  #scale_x_discrete(
  #breaks = c("Primordia", "Differentiation", "Fundamental", "Late fruiting", "Late non fruiting"),
  #labels = c("Primordia", "Differentiation", "Fundamental", "Late fruiting", "Late non fruiting"), 
  #drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class > 2%) \n") 
plot(bar_fungi)

#Sordariomycetes only
tax_table(biom_ITS_uparse_in)
sord_only = subset_taxa(biom_ITS_uparse_in, Class=="Sordariomycetes")
sord_only
tax_table(sord_only)
write.csv(tax_table(sord_only), "sord_only.csv")

# taxonomy table was exported and the unidentified genera that could be identified using NCBI BLAST were corrected
tax_ITS_in_correct = read.csv("sord_only_corrected.csv", header=T, row.names =1)
tax_table(sord_only) <- tax_table(as.matrix(tax_ITS_in_correct))
tax_table(sord_only) 
biom_ITS_uparse_sord_genus <- sord_only %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by genus

bar_fungi_sord = ggplot(biom_ITS_uparse_sord_genus, aes(x = Label, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits=c("Primordia 1", "Primordia 2", "Primordia 3", "Differentiation 1", "Differentiation 2","Differentiation 3", "Fundaments 1","Fundaments 2", "Fundaments 3", "Late fruiting 1", "Late fruiting 2", "Late fruiting 3", "Late non fruiting 1", "Late non fruiting 2", "Late non fruiting 3"))+
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 2%) \n") 
bar_fungi_sord


#combining stacked barplots (FIGURE 1)---------------
ggarrange(bar_prok,
          bar_fungi,
          bar_fungi_sord,
          labels = c("A","B","C"),
          label.x = 0, label.y = 1, hjust = -1.0, vjust = 1.5,
          widths = c(1, 1),
          align = "v",
          ncol = 3, nrow = 1 )





# alpha diversity ------------
library("ggpubr")
library("gridExtra")
library("grid")
library("cowplot")

# prokaryotes
sample_data(biom_16s_uparse_in)
levels(sample_data(biom_16s_uparse_in)$Stage)[levels
                                              (sample_data(biom_16s_uparse_in)$Stage)=="Late_fruiting"] <- "Late fruiting"

levels(sample_data(biom_16s_uparse_in)$Stage)[levels
                                              (sample_data(biom_16s_uparse_in)$Stage)=="Late_non_fruiting"] <- "Late non fruiting"


sample_data(biom_16s_uparse_in)$Stage <- factor(sample_data(biom_16s_uparse_in)$Stage,
                                                levels=c("Primordia","Fundaments","Differentiation","Late fruiting","Late non fruiting"))

sample_data(biom_16s_uparse_in)
alpha_16s_in = plot_richness(biom_16s_uparse_in, x="Stage", color = "Stage",
                             measures = c("Shannon", "Observed")) +
  theme_pubclean() +
  geom_point(size = 0, shape = 16) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  expand_limits(x = 0, y = 0) +
  labs(title="", x="Stage", y = "") +
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF"))+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())

alpha_16s_in
alpha_16s_in$layers
alpha_16s_in$layers <- alpha_16s_in$layers[-1]
alpha_16s_in$layers <- alpha_16s_in$layers[-1]
alpha_16s_in

# rename facets on the plot 
levels(alpha_16s_in$data$variable)
levels(alpha_16s_in$data$variable)[levels
                                   (alpha_fungi_in$data$variable)=="Observed"] <- "Richness"
alpha_16s_in


# fungi
sample_data(biom_ITS_uparse_in)
levels(sample_data(biom_ITS_uparse_in)$Stage)[levels
                                              (sample_data(biom_ITS_uparse_in)$Stage)=="Late_fruiting"] <- "Late fruiting"

levels(sample_data(biom_ITS_uparse_in)$Stage)[levels
                                              (sample_data(biom_ITS_uparse_in)$Stage)=="Late_non_fruiting"] <- "Late non fruiting"


sample_data(biom_ITS_uparse_in)$Stage <- factor(sample_data(biom_ITS_uparse_in)$Stage,
                                                levels=c("Primordia","Fundaments","Differentiation","Late fruiting","Late non fruiting"))


alpha_fungi_in = plot_richness(biom_ITS_uparse_in, x="Stage", color = "Stage",
                               measures = c("Shannon", "Observed")) +
  theme_pubclean() +
  geom_point(size = 0, shape = 16) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  expand_limits(x = 0, y = 0) +
  labs(title="", x="Stage", y = "") +
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())

alpha_fungi_in
alpha_fungi_in$layers
alpha_fungi_in$layers <- alpha_fungi_in$layers[-1]
alpha_fungi_in$layers <- alpha_fungi_in$layers[-1]
alpha_fungi_in

# rename facets on the plot 
levels(alpha_fungi_in$data$variable)
levels(alpha_fungi_in$data$variable)[levels
                                     (alpha_fungi_in$data$variable)=="Observed"] <- "Richness"
alpha_fungi_in

#create combined plot (FIGURE 2)----------------
library("ggpubr")

ggarrange(alpha_16s_in,
          alpha_fungi_in,
          labels = c("A","B"),
          widths = c(1, 1),
          align = "v",
          ncol = 2, nrow = 1 )

###creating Venn Diagrams-----------------

library("limma")

#source("Rscript_1_venn_diagram.R")

# my.venn_Gian.R modified form limma package
# function for venn
my.venn_Gian <- function (object, include = "both", names = NULL, mar = rep(1,
                                                                            4), cex = c(1.5, 1, 0.7), lwd = 1, circle.col = NULL, counts.col = NULL,
                          show.include = NULL, ...)
{
  include <- as.character(include)
  LenInc <- min(length(include), 2)
  if (is(object, "VennCounts")) {
    include <- include[1]
    LenInc <- 1
  }
  else {
    if (LenInc > 1)
      z2 <- vennCounts(object, include = include[2])[,
                                                     "Counts"]
    object <- vennCounts(object, include = include[1])
  }
  z <- object[, "Counts"]
  nsets <- ncol(object) - 1
  if (nsets > 5)
    stop("Can't plot Venn diagram for more than 5 sets")
  VennZone <- object[, 1:nsets, drop = FALSE]
  VennZone <- apply(VennZone, 1, function(x) paste(x, sep = "",
                                                   collapse = ""))
  names(z) <- VennZone
  if (length(include) == 2)
    names(z2) <- VennZone
  if (is.null(names))
    names <- colnames(object)[1:nsets]
  FILL.COL <- TRUE
  if (is.null(circle.col)) {
    circle.col <- par("col")
    FILL.COL <- FALSE
  }
  if (length(circle.col) < nsets)
    circle.col <- rep(circle.col, length.out = nsets)
  if (is.null(counts.col))
    counts.col <- par("col")
  if (length(counts.col) < LenInc)
    counts.col <- rep(counts.col, length.out = LenInc)
  if (is.null(show.include))
    show.include <- as.logical(LenInc - 1)
  old.par <- par()$mar
  on.exit(par(mar = old.par))
  par(mar = mar)
  if (nsets <= 3) {
    plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4,
                                                             4), xlab = "", ylab = "", axes = FALSE, ...)
    theta <- 2 * pi * (0:360)/360
    xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
    ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
    r <- 1.5
    xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2,
                                                 0))
    ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(2.4, 2.4,
                                               -3))
    for (circle in 1:nsets) {
      if (!FILL.COL)
        lines(xcentres[circle] + r * cos(theta), ycentres[circle] +
                r * sin(theta), lwd = lwd, col = circle.col[circle])
      if (FILL.COL) {
        RGB <- col2rgb(circle.col[circle])/255
        ALPHA <- 0.06
        RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                       alpha = ALPHA)
        polygon(xcentres[circle] + r * cos(theta), ycentres[circle] +
                  r * sin(theta), border = circle.col[circle],
                lwd = lwd, col = RGB.ALP)
      }
      text(xtext[circle], ytext[circle], names[circle],
           cex = cex)
    }
    #switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5,
    #    3, 2.5), rect(-3, -3.5, 3, 3.3))
    showCounts <- switch(nsets, function(counts, cex, adj,
                                         col, leg) {
      text(2.3, -2.1, sprintf("Not in any =  %i", counts[1]), cex = cex, col = col,
           adj = adj)
      text(0, 0, counts[2], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.3, -2.1, leg, cex = cex,
                             col = col, adj = adj)
    }, function(counts, cex, adj, col, leg) {
      text(2.3, -2.1, sprintf("Not in any = %i", counts[1]), cex = cex, col = col,
           adj = adj)
      text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
      text(-1.5, 0.1, counts[3], cex = cex, col = col,
           adj = adj)
      text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.3, -2.1, leg, cex = cex,
                             col = col, adj = adj)
    }, function(counts, cex, adj, col, leg) {
      text(2.5, -3, sprintf("Not in any = %i", counts[1]), cex = cex, col = col, adj = adj)
      text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
      text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
      text(0.75, -0.35, counts[4], cex = cex, col = col,
           adj = adj)
      text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
      text(-0.75, -0.35, counts[6], cex = cex, col = col,
           adj = adj)
      text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
      text(0, 0, counts[8], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.5, -3, leg, cex = cex,
                             col = col, adj = adj)
    })
    if (LenInc == 1)
      adj <- c(0.5, 0.5)
    else adj <- c(0.5, 0)
    print(z)
    showCounts(counts = z, cex = cex[1], adj = adj, col = counts.col[1],
               leg = include[1])
    if (LenInc == 2)
      showCounts(counts = z2, cex = cex[1], adj = c(0.5,
                                                    1), col = counts.col[2], leg = include[2])
    return(invisible())
  }
  plot(c(-20, 420), c(-20, 420), type = "n", axes = FALSE,
       ylab = "", xlab = "", ...)
  relocate_elp <- function(e, alpha, x, y) {
    phi <- (alpha/180) * pi
    xr <- e[, 1] * cos(phi) + e[, 2] * sin(phi)
    yr <- -e[, 1] * sin(phi) + e[, 2] * cos(phi)
    xr <- x + xr
    yr <- y + yr
    cbind(xr, yr)
  }
  if (4 == nsets) {
    #rect(-20, -20, 420, 400)
    elps <- cbind(162 * cos(seq(0, 2 * pi, len = 1000)),
                  108 * sin(seq(0, 2 * pi, len = 1000)))
    if (!FILL.COL) {
      polygon(relocate_elp(elps, 45, 130, 170), border = circle.col[1],
              lwd = lwd)
      polygon(relocate_elp(elps, 45, 200, 200), border = circle.col[2],
              lwd = lwd)
      polygon(relocate_elp(elps, 135, 200, 200), border = circle.col[3],
              lwd = lwd)
      polygon(relocate_elp(elps, 135, 270, 170), border = circle.col[4],
              lwd = lwd)
    }
    if (FILL.COL) {
      RGB <- col2rgb(circle.col)/255
      ALPHA <- 0.06
      RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                      alpha = ALPHA)
      RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2],
                      alpha = ALPHA)
      RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3],
                      alpha = ALPHA)
      RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4],
                      alpha = ALPHA)
      polygon(relocate_elp(elps, 45, 130, 170), border = circle.col[1],
              lwd = lwd, col = RGB.ALP1)
      polygon(relocate_elp(elps, 45, 200, 200), border = circle.col[2],
              lwd = lwd, col = RGB.ALP2)
      polygon(relocate_elp(elps, 135, 200, 200), border = circle.col[3],
              lwd = lwd, col = RGB.ALP3)
      polygon(relocate_elp(elps, 135, 270, 170), border = circle.col[4],
              lwd = lwd, col = RGB.ALP4)
    }
    text(35, 315, names[1], cex = cex[1])
    text(138, 350, names[2], cex = cex[1])
    text(262, 347, names[3], cex = cex[1])
    text(365, 315, names[4], cex = cex[1])
    text(35, 250, z["1000"], cex = cex[2], col = counts.col[1],
    )
    text(140, 315, z["0100"], cex = cex[2], col = counts.col[1])
    text(260, 315, z["0010"], cex = cex[2], col = counts.col[1])
    text(365, 250, z["0001"], cex = cex[2], col = counts.col[1])
    text(90, 282, z["1100"], cex = cex[3], col = counts.col[1])
    text(95, 110, z["1010"], cex = cex[2], col = counts.col[1])
    text(200, 52, z["1001"], cex = cex[3], col = counts.col[1])
    text(200, 292, z["0110"], cex = cex[2], col = counts.col[1])
    text(300, 110, z["0101"], cex = cex[2], col = counts.col[1])
    text(310, 282, z["0011"], cex = cex[3], col = counts.col[1])
    text(130, 230, z["1110"], cex = cex[2], col = counts.col[1])
    text(245, 81, z["1101"], cex = cex[3], col = counts.col[1])
    text(155, 81, z["1011"], cex = cex[3], col = counts.col[1])
    text(270, 230, z["0111"], cex = cex[2], col = counts.col[1])
    text(200, 152, z["1111"], cex = cex[2], col = counts.col[1])
    text(400, 15, sprintf("Not in any = %i", z["0000"]), cex = cex[1], col = counts.col[1])
    if (length(include) == 2) {
      text(35, 238, z2["1000"], cex = cex[2], col = counts.col[2])
      text(140, 304, z2["0100"], cex = cex[2], col = counts.col[2])
      text(260, 304, z2["0010"], cex = cex[2], col = counts.col[2])
      text(365, 238, z2["0001"], cex = cex[2], col = counts.col[2])
      text(90, 274, z2["1100"], cex = cex[3], col = counts.col[2])
      text(95, 100, z2["1010"], cex = cex[2], col = counts.col[2])
      text(200, 43, z2["1001"], cex = cex[3], col = counts.col[2])
      text(200, 280, z2["0110"], cex = cex[2], col = counts.col[2])
      text(300, 100, z2["0101"], cex = cex[2], col = counts.col[2])
      text(310, 274, z2["0011"], cex = cex[3], col = counts.col[2])
      text(130, 219, z2["1110"], cex = cex[2], col = counts.col[2])
      text(245, 71, z2["1101"], cex = cex[3], col = counts.col[2])
      text(155, 72, z2["1011"], cex = cex[3], col = counts.col[2])
      text(270, 219, z2["0111"], cex = cex[2], col = counts.col[2])
      text(200, 140, z2["1111"], cex = cex[2], col = counts.col[2])
      text(400, -2, sprintf("Not in any = %i", z2["0000"]), cex = cex[1], col = counts.col[2])
      if (show.include) {
        text(10, 15, include[1], cex = cex[1], col = counts.col[1])
        text(10, -2, include[2], cex = cex[1], col = counts.col[2])
      }
    }
    return(invisible())
  }
  #rect(-20, -30, 430, 430)
  elps <- cbind(150 * cos(seq(0, 2 * pi, len = 1000)), 60 *
                  sin(seq(0, 2 * pi, len = 1000)))
  if (!FILL.COL) {
    polygon(relocate_elp(elps, 90, 200, 250), border = circle.col[1],
            lwd = lwd)
    polygon(relocate_elp(elps, 162, 250, 220), border = circle.col[2],
            lwd = lwd)
    polygon(relocate_elp(elps, 234, 250, 150), border = circle.col[3],
            lwd = lwd)
    polygon(relocate_elp(elps, 306, 180, 125), border = circle.col[4],
            lwd = lwd)
    polygon(relocate_elp(elps, 378, 145, 200), border = circle.col[5],
            lwd = lwd)
  }
  if (FILL.COL) {
    RGB <- col2rgb(circle.col)/255
    ALPHA <- 0.06
    RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
    RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2], alpha = ALPHA)
    RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3], alpha = ALPHA)
    RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4], alpha = ALPHA)
    RGB.ALP5 <- rgb(RGB[1, 5], RGB[2, 5], RGB[3, 5], alpha = ALPHA)
    polygon(relocate_elp(elps, 90, 200, 250), border = circle.col[1],
            lwd = lwd, col = RGB.ALP1)
    polygon(relocate_elp(elps, 162, 250, 220), border = circle.col[2],
            lwd = lwd, col = RGB.ALP2)
    polygon(relocate_elp(elps, 234, 250, 150), border = circle.col[3],
            lwd = lwd, col = RGB.ALP3)
    polygon(relocate_elp(elps, 306, 180, 125), border = circle.col[4],
            lwd = lwd, col = RGB.ALP4)
    polygon(relocate_elp(elps, 378, 145, 200), border = circle.col[5],
            lwd = lwd, col = RGB.ALP5)
  }
  text(50, 285, names[1], cex = cex[1])
  text(200, 415, names[2], cex = cex[1])
  text(350, 305, names[3], cex = cex[1])
  text(350, 20, names[4], cex = cex[1])
  text(100, -10, names[5], cex = cex[1])
  text(61, 231, z["10000"], cex = cex[2], col = counts.col[1])
  text(200, 332, z["01000"], cex = cex[2], col = counts.col[1])
  text(321, 248, z["00100"], cex = cex[2], col = counts.col[1])
  text(290, 84, z["00010"], cex = cex[2], col = counts.col[1])
  text(132, 72, z["00001"], cex = cex[2], col = counts.col[1])
  text(146, 253, z["11000"], cex = cex[3], col = counts.col[1])
  text(123, 191, z["10100"], cex = cex[3], col = counts.col[1])
  text(275, 155, z["10010"], cex = cex[3], col = counts.col[1])
  text(137, 149, z["10001"], cex = cex[3], col = counts.col[1])
  text(243, 271, z["01100"], cex = cex[3], col = counts.col[1])
  text(175, 270, z["01010"], cex = cex[3], col = counts.col[1])
  text(187, 120, z["01001"], cex = cex[3], col = counts.col[1])
  text(286, 193, z["00110"], cex = cex[3], col = counts.col[1])
  text(267, 238, z["00101"], cex = cex[3], col = counts.col[1])
  text(228, 108, z["00011"], cex = cex[3], col = counts.col[1])
  text(148, 213, z["11100"], cex = cex[3], col = counts.col[1])
  text(159, 255, z["11010"], cex = cex[3], col = counts.col[1])
  text(171, 144, z["11001"], cex = cex[3], col = counts.col[1])
  text(281, 178, z["10110"], cex = cex[3], col = counts.col[1])
  text(143, 166, z["10101"], cex = cex[3], col = counts.col[1])
  text(252, 148, z["10011"], cex = cex[3], col = counts.col[1])
  text(205, 258, z["01110"], cex = cex[3], col = counts.col[1])
  text(254, 248, z["01101"], cex = cex[3], col = counts.col[1])
  text(211, 121, z["01011"], cex = cex[3], col = counts.col[1])
  text(267, 214, z["00111"], cex = cex[3], col = counts.col[1])
  text(170, 234, z["11110"], cex = cex[3], col = counts.col[1])
  text(158, 172, z["11101"], cex = cex[3], col = counts.col[1])
  text(212, 142, z["11011"], cex = cex[3], col = counts.col[1])
  text(263, 183, z["10111"], cex = cex[3], col = counts.col[1])
  text(239, 235, z["01111"], cex = cex[3], col = counts.col[1])
  text(204, 193, z["11111"], cex = cex[2], col = counts.col[1])
  text(400, 7, sprintf("Not in any = %i", z["00000"]), cex = cex[1], col = counts.col[1])
  if (length(include) == 2) {
    text(61, 220, z2["10000"], cex = cex[2], col = counts.col[2])
    text(200, 321, z2["01000"], cex = cex[2], col = counts.col[2])
    text(321, 237, z2["00100"], cex = cex[2], col = counts.col[2])
    text(290, 73, z2["00010"], cex = cex[2], col = counts.col[2])
    text(132, 61, z2["00001"], cex = cex[2], col = counts.col[2])
    text(146, 244, z2["11000"], cex = cex[3], col = counts.col[2])
    text(123, 180, z2["10100"], cex = cex[3], col = counts.col[2])
    text(275, 144, z2["10010"], cex = cex[3], col = counts.col[2])
    text(137, 143, z2["10001"], cex = cex[3], col = counts.col[2])
    text(243, 260, z2["01100"], cex = cex[3], col = counts.col[2])
    text(175, 259, z2["01010"], cex = cex[3], col = counts.col[2])
    text(187, 110, z2["01001"], cex = cex[3], col = counts.col[2])
    text(286, 186, z2["00110"], cex = cex[3], col = counts.col[2])
    text(267, 230, z2["00101"], cex = cex[3], col = counts.col[2])
    text(228, 97, z2["00011"], cex = cex[3], col = counts.col[2])
    text(148, 203, z2["11100"], cex = cex[3], col = counts.col[2])
    text(159, 249, z2["11010"], cex = cex[3], col = counts.col[2])
    text(171, 137, z2["11001"], cex = cex[3], col = counts.col[2])
    text(281, 171, z2["10110"], cex = cex[3], col = counts.col[2])
    text(143, 155, z2["10101"], cex = cex[3], col = counts.col[2])
    text(252, 137, z2["10011"], cex = cex[3], col = counts.col[2])
    text(205, 247, z2["01110"], cex = cex[3], col = counts.col[2])
    text(254, 242, z2["01101"], cex = cex[3], col = counts.col[2])
    text(211, 112, z2["01011"], cex = cex[3], col = counts.col[2])
    text(267, 207, z2["00111"], cex = cex[3], col = counts.col[2])
    text(170, 223, z2["11110"], cex = cex[3], col = counts.col[2])
    text(158, 162, z2["11101"], cex = cex[3], col = counts.col[2])
    text(212, 133, z2["11011"], cex = cex[3], col = counts.col[2])
    text(263, 172, z2["10111"], cex = cex[3], col = counts.col[2])
    text(239, 228, z2["01111"], cex = cex[3], col = counts.col[2])
    text(204, 182, z2["11111"], cex = cex[2], col = counts.col[2])
    text(400, -10, sprintf("Not in any = %i", z2["00000"]), cex = cex[1], col = counts.col[2])
    if (show.include) {
      text(10, 7, include[1], cex = cex[1], col = counts.col[1])
      text(10, -10, include[2], cex = cex[1], col = counts.col[2])
    }
  }
  invisible()
}




# fungi venn diagram

biom_ITS_uparse_in_St = merge_samples(biom_ITS_uparse_in, "Stage")
otu_fungi_in_St <- as.data.frame(t(otu_table(biom_ITS_uparse_in_St)))
venn_counts_otu_fungi_in_St <- vennCounts(otu_fungi_in_St, include="both")
venn_counts_otu_fungi_in_St


biom_16s_uparse_in_st = merge_samples(biom_16s_uparse_in, "Stage")
otu_prokaryote_in_St <- as.data.frame(t(otu_table(biom_16s_uparse_in_st)))
venn_counts_otu_prokaryote_in_St <- vennCounts(otu_prokaryote_in_St, include="both")
venn_counts_otu_prokaryote_in_St


# *** Creating Venn Diagrams (FIGURE 4) --------------------------------

layout(matrix(1:2, ncol=2))


prok_venn = my.venn_Gian(venn_counts_otu_prokaryote_in_St,
                         cex=c(0.9), 
                         circle.col = "blue",
                         mar = c(1,1,1,1),
                         lwd = 2, main="A")
fung_venn = my.venn_Gian(venn_counts_otu_fungi_in_St,
                         cex=c(0.9), 
                         circle.col = "black",
                         mar = c(1,1,1,1),
                         lwd = 2, main="B")



dev.off()

### Normalizing Data with Metagenomeseq------------------
source("http://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
library("metagenomeSeq")


# fitting into a Gaussian Model using metagenomeSeq 16S
otu_table(biom_16s_uparse_in)
biom_16s_uparse_in_norm = phyloseq_to_metagenomeSeq(biom_16s_uparse_in)
otu_table(biom_16s_uparse_in_norm)
p_biom<-cumNormStatFast(biom_16s_uparse_in_norm)
biom_quant<-cumNorm(biom_16s_uparse_in_norm, p=p_biom)
biom_quant
normFactors(biom_quant)
biom_16s_uparse_in_norm<-MRcounts(biom_quant, norm=T)
head(biom_16s_uparse_in_norm)
biom_16s_uparse_in_norm
#create physeq object with normalized otu table
otu_table(biom_16s_uparse_in) <- otu_table(biom_16s_uparse_in_norm, taxa_are_rows = TRUE)


# fitting into a Gaussian Model using metagenomeSeq ITS
otu_table(biom_ITS_uparse_in)
biom_ITS_uparse_in_norm = phyloseq_to_metagenomeSeq(biom_ITS_uparse_in)
otu_table(biom_ITS_uparse_in_norm)
p_biom<-cumNormStatFast(biom_ITS_uparse_in_norm)
biom_quant<-cumNorm(biom_ITS_uparse_in_norm, p=p_biom)
biom_quant
normFactors(biom_quant)
biom_ITS_uparse_in_norm<-MRcounts(biom_quant, norm=T)
head(biom_ITS_uparse_in_norm)
biom_ITS_uparse_in_norm
#create physeq object with normalized otu table
otu_table(biom_ITS_uparse_in) <- otu_table(biom_ITS_uparse_in_norm, taxa_are_rows = TRUE)
otu_table(biom_ITS_uparse_in)

### Beta Diversity - NMDS plots-----------------

# prokaryotes
nmds_16s_in = ordinate(biom_16s_uparse_in, method ="NMDS", distance="bray", try=200)
nmds_16s_in

p_nmds_16s_in = plot_ordination(biom_16s_uparse_in, nmds_16s_in, color="Stage", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, bg= "#0000FF") + # ,aes(shape=Age))
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF")) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  theme(legend.position="right")

p_nmds_16s_in


# fungi

p_nmds_ITS_in = plot_ordination(biom_ITS_uparse_in, nmds_ITS_in, color="Stage", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, bg="black") + # ,aes(shape=Age))
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000")) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #theme(legend.title=element_blank()) + 
  theme(legend.position="right")
plot(p_nmds_ITS_in)

### Combining plots (FIGURE 3)----------
ggarrange(p_nmds_16s_in,
          p_nmds_ITS_in,
          labels = c("A","B"),
          widths = c(2, 2),
          align = "v",
          ncol = 2, nrow = 1 )


### PERMANOVA (TABLE 1)-------

# creating vegan objects for use in permanova

otu_fungi_in <- as.data.frame(otu_table(biom_ITS_uparse_in))
taxa_fungi_in <- as.data.frame(as.matrix(tax_table(biom_ITS_uparse_in)))
metadata_fungi_in <- as.data.frame(as.matrix(sample_data(biom_ITS_uparse_in)))
dim(otu_fungi_in)



otu_prokaryote_in <- as.data.frame(otu_table(biom_16s_uparse_in))
taxa_prokaryote_in <- as.data.frame(as.matrix(tax_table(biom_16s_uparse_in)))
metadata_prokaryote_in <- as.data.frame(as.matrix(sample_data(biom_16s_uparse_in)))
dim(otu_prokaryote_in)
otu_prokaryote_in

# permanova prokaryote
metadata_prokaryote_in
adonis(t(otu_prokaryote_in) ~ Stage, data=metadata_prokaryote_in, permutations=9999)

vegan::vegdist(t(otu_prokaryote_in), method="bray") -> dist_otu_prokaryote_in
# betadisper prokaryote
permdisp_otu_prokaryote_in_S <- betadisper(dist_otu_prokaryote_in, metadata_prokaryote_in$Stage)
anova(permdisp_otu_prokaryote_in_S)

#PERMANOVA fungi
metadata_fungi_in
adonis(t(otu_fungi_in) ~ Stage, data=metadata_fungi_in, permutations=9999)

vegan::vegdist(t(otu_fungi_in), method="bray") -> dist_otu_fungi_in

# betadisper fungi
permdisp_otu_fungi_in_S <- betadisper(dist_otu_fungi_in, metadata_fungi_in$Stage)
anova(permdisp_otu_fungi_in_S)

