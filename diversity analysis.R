library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names



#dummy data for practise
data("GlobalPatterns")

GlobalPatterns
GlobalPatterns |> 
  otu_table() |> 
  view()

# transforming to phyloseq object -----------------------------------------

# importing the required tables
otu_mat<- read_excel("Anemia Assignment Data OTU.xlsx", sheet = "Anemia Assignment Data OTU")
tax_mat<- read_excel("Anemia Assignment Data OTU.xlsx", sheet = "taxonomy")
samples_df <- read_excel("Anemia Assignment Data OTU.xlsx", sheet = "metadata")

#creating rownames for all datafiles of phyloseq object
otu_mat <- otu_mat |> 
  column_to_rownames(var = 'OTU')
tax_mat <- tax_mat |> 
  column_to_rownames(var = 'OTU')
samples_df <- samples_df |> 
  column_to_rownames(var = 'Sample_ID')

#tax and otu table needed to be in matrix format
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#transforming to phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

carbom <- phyloseq(OTU, TAX, samples)
carbom

# visualising data
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

#alpha diversity based on anemia


plot_richness(carbom, x="ANEMIA..YES.NO.", measures=c("Observed", "Shannon","Chao1",'Simpson')) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

#p value using wilklnsin test
rich = estimate_richness(carbom, measures = c("Observed", "Shannon","Chao1",'Simpson'))
wilcox.observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(carbom)$`ANEMIA (YES/NO)`, 
                                        p.adjust.method = "BH")

tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.observed


wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(carbom)$`ANEMIA (YES/NO)`, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.shannon


wilcox.Chao1 <- pairwise.wilcox.test(rich$Chao1, 
                                       sample_data(carbom)$`ANEMIA (YES/NO)`, 
                                       p.adjust.method = "BH")
tab.Chao1 <- wilcox.Chao1$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.Chao1


wilcox.Simpson <- pairwise.wilcox.test(rich$Simpson, 
                                     sample_data(carbom)$`ANEMIA (YES/NO)`, 
                                     p.adjust.method = "BH")
tab.Simpson <- wilcox.Simpson$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.Simpson


#beta analysis

dist = phyloseq::distance(carbom, method="bray")
ordination = ordinate(carbom, method="PCoA", distance=dist)
plot_ordination(carbom, ordination, color="ANEMIA..YES.NO.") + 
  theme_classic() +
  theme(strip.background = element_blank()) 



pcoa <- ordinate(carbom, method = "PCoA", distance = dist)

plot_ordination(carbom, pcoa) +
  ordisurf(pcoa, NLR, add = TRUE)
