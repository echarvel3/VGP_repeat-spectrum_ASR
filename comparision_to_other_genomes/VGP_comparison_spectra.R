
library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)

setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/comparision_to_other_genomes/")
source("./subset_phylogeny.R")

# functions 
split_row <- function(x){
  z = paste(collapse = "_", strsplit(x, split = "_")[[1]][1:2])
}

# reading spectrum data for existing genomes
low_qual_spec_params = merge(
  fread("./respect_full-existing_spectra.txt"),
  fread("./respect_full-existing_parameters.txt")
)
low_qual_spec = merge(low_qual_spec_params, existing_merged[,c("ScientificName", "accession", "sample")])


high_qual_spec_params = merge(
  fread("../full_VGP_analyses/respect_full-VGP_spectra.txt"), 
  fread("../full_VGP_analyses/respect_full-VGP_parameters.txt")
)
high_qual_spec = merge(high_qual_spec_params, vgp_merged[,c("ScientificName", "sample", "Extended lineage")])



# formatting vgp data
high_qual_spec[,2:51] = high_qual_spec[,2:51]/1000
high_qual_spec$genome_length = high_qual_spec$genome_length/1000
x = melt(data = high_qual_spec, id.vars = c("ScientificName", "Extended lineage", "genome_length"), measure.vars = c(2:51))


# formatting existing data
low_qual_spec[,2:51] = low_qual_spec[,2:51]/1000
low_qual_spec$genome_length = as.numeric(low_qual_spec$genome_length)/1000
y = melt(data = low_qual_spec, id.vars = c("ScientificName", "genome_length"), measure.vars = c(2:51))

z = merge(x=x, y=y, by = c("ScientificName", "variable"))
z$i = as.numeric(gsub("r", "", z$variable))
z = z[!ScientificName %in%  c("Oenanthe_melanoleuca", "Rhinophrynus_dorsalis") & ScientificName %in% tree$tip.label]
z$ratio = z$value.x/z$value.y
z$ratio

head(z)
  library(dplyr)
library(scales)
arrow_plot = z[!is.na(`Extended lineage`)] %>%
  group_by(`Extended lineage`,i) %>%
  summarize(newm = as.numeric(mean(value.y*i)),oldm = as.numeric(mean(value.x*i))) %>%
  #summarize(newm = mean(value.y*i/genome_length.y),oldm = mean(value.x*i/genome_length.x)) %>%
  
ggplot(
                 aes(
                   x = i,
                   #y = sin_y((value.x - value.y)/value.y),
                   yend =  as.numeric(newm*1000),
                   y = as.numeric(oldm*1000),
                   color = `Extended lineage`)) +
  geom_segment(size = 0.65, alpha = 1,position = position_dodge(width=0.5),
               arrow = arrow(ends = "first",length =  unit(4,"pt"))) +
  # facet_wrap(.~ifelse(
  #   `Extended lineage` %in% c("Lepidosauria", "Turtles", "Crocodiles", "Birds", "Mammals"), 
  #   "Amniotes", 
  #   ifelse(`Extended lineage` %in% c("Lobe-finned fishes", "Ray-finned fishes", "Cartilaginous fishes", "Cyclostomes"), "Fish", "Amphibians")), 
  #   nrow = 1, scales = "free_y") +
  #facet_wrap(.~i<3, scales = "free") +
  theme_classic() + 
  scale_color_manual(name=NULL, values = colorRampPalette(RColorBrewer::brewer.pal(n=8, name = "Accent"))(10)) +
   scale_color_manual(name=NULL,
    values = c("Mammals" = "#E69F00",  # Okabe–Ito Orange
                                "Birds" = "#00796B",    # Teal 700
                                "Crocodiles" = "#009688",# Teal 500
                                "Turtles" = "#4DB6AC",  # Teal 300
                                "Lepidosauria" = "#80CBC4", # Teal 200
                                "Amphibians" = "#984EA3",   # Purple
                                "Lobe-finned fishes" = "#A6761D", # Brown
                                "Ray-finned fishes" = "#56B4E9",  # Okabe–Ito Sky Blue
                                "Cartilaginous fishes" = "#0072B2", # Okabe–Ito Blue
                                "Cyclostomes" = "#CC79A7", # Okabe–Ito Magenta
                                "Ancestral" = "#999999"  # Neutral Gray
  )) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        title = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
        legend.position = c(0.9345, 0.79),
        legend.box.background = element_rect(color="black", size=1),
        #legend.justification = c("right", "top"),
        legend.box.just = c("right", "top"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(name="k-mer count",
                     #labels=percent,
                     trans="log10")+
  scale_x_continuous(breaks=c(1,10,20,30,40,50),name="31-mer frequency")
#arrow_plot
ggsave("../figures/arrow_plot_kmer_count.pdf", plot = arrow_plot, width = 19, height = 7)

#  Removed kmer frequency plot
#  z[!is.na(`Extended lineage`)] %>%
#   arrange(`ScientificName`,-i) %>%
#   group_by(`ScientificName`) %>% #print(n=70) 
#   mutate(newm = as.numeric(cumsum(value.y*i)),
#             oldm = as.numeric(cumsum(value.x*i))) %>% #print(n=70) 
#    group_by(`Extended lineage`,i) %>%
#    summarize(newmcc = mean(newm),oldmcc = mean(oldm)) %>% #print(n=70) 
# ggplot(
#   aes(
#     x = i,
#     #y = sin_y((value.x - value.y)/value.y),
#     yend =  as.numeric(newmcc*1000),
#     y = as.numeric(oldmcc*1000),
#     color = `Extended lineage`)) +
#   geom_segment(size = 0.65, alpha = 1,position = position_dodge(width=0.5),
#                arrow = arrow(ends = "first",length =  unit(4,"pt"))) +
#    facet_wrap(.~ifelse(
#      `Extended lineage` %in% c("Ray-finned fishes", "Cyclostomes", "Birds","Lepidosauria"), 
#      "Small", 
#      ifelse(`Extended lineage` %in% c("Lobe-finned fishes"), "Large", "Medium")), 
#      # "Ray-finned fishes", "Cartilaginous fishes", "Cyclostomes"
#      nrow = 1, scales = "free_y") +
#   #facet_wrap(.~i<3, scales = "free") +
#   theme_classic() + 
#   scale_color_manual(name=NULL, values = colorRampPalette(RColorBrewer::brewer.pal(n=8, name = "Accent"))(10)) +
#   scale_color_manual(name=NULL,
#                      values = c("Mammals" = "#E69F00",  # Okabe–Ito Orange
#                                 "Birds" = "#00796B",    # Teal 700
#                                 "Crocodiles" = "#009688",# Teal 500
#                                 "Turtles" = "#4DB6AC",  # Teal 300
#                                 "Lepidosauria" = "#80CBC4", # Teal 200
#                                 "Amphibians" = "#984EA3",   # Purple
#                                 "Lobe-finned fishes" = "#A6761D", # Brown
#                                 "Ray-finned fishes" = "#56B4E9",  # Okabe–Ito Sky Blue
#                                 "Cartilaginous fishes" = "#0072B2", # Okabe–Ito Blue
#                                 "Cyclostomes" = "#CC79A7", # Okabe–Ito Magenta
#                                 "Ancestral" = "#999999"  # Neutral Gray
#                      )) +
#   theme(axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size=15),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         title = element_text(size=15),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
#         legend.position =  "bottom",
#         legend.box.background = element_rect(color="black", size=1),
#         #legend.justification = c("right", "top"),
#         legend.box.just = c("right", "top"),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=20),
#         legend.key.spacing.x = unit(2.0, "cm"),
#         strip.text = element_text(size = 17),
#         panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#         plot.margin = margin(10, 10, 10, 10)) +
#   scale_y_continuous(name="k-mer count",
#                      #labels=percent,
#                      trans="log10")+
#   scale_x_continuous(breaks=c(1,10,20,30,40,50),name="31-mer frequency")
