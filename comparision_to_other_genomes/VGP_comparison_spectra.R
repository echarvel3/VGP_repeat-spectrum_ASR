
library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)

setwd("~/Desktop/VGP_analyses/")

split_row <- function(x){
  z = paste(collapse = "_", strsplit(x, split = "_")[[1]][1:2])
}

low_qual_spec = fread("./low_qual_genomes/vertebrate_querySpeciesRep.tsv")

low_qual_spec$sample = low_qual_spec$assembly_id
low_qual_spec = low_qual_spec[,!c("assembly_id")]

low_qual_spec[scientific_name == "Gorilla gorilla gorilla", scientific_name := "Gorilla gorilla"]
low_qual_spec[scientific_name == "Phoenicopterus ruber ruber", scientific_name := "Phoenicopterus ruber"]
low_qual_spec[scientific_name == "Canis lupus", scientific_name := "XXXX"]
low_qual_spec[scientific_name == "Alca torda", scientific_name := "Alca Torda"]
low_qual_spec[scientific_name == "Mustela nivalis", scientific_name := "Mustela nivalis vulgaris"]
low_qual_spec[scientific_name == "Monodon monoceros", scientific_name := "Monodon monocero"]
low_qual_spec[scientific_name == "Lonchura striata", scientific_name := "Lonchura striata domestica"]
low_qual_spec[scientific_name == "Guaruba guarouba", scientific_name := "Guaruba guaruba"]
low_qual_spec[scientific_name == "Salvelinus alpinus oquassa", scientific_name := "Salvelinus alpinus"]
low_qual_spec[scientific_name == "Rhamphochromis sp. 'chilingali'", scientific_name := "Rhamphochromis sp. chilingali"]
low_qual_spec[scientific_name == "Lycocorax pyrrhopterus obiensis", scientific_name := "Lycocorax pyrrhopterus"]

low_qual_spec_params = fread("./VGP_repeat-spectrum_ASR/comparision_to_other_genomes/respect_full-existing_spectra.txt")
low_qual_spec_params$sample = sapply(X = low_qual_spec_params$sample, FUN = split_row)
low_qual_spec_params_new = merge(low_qual_spec_params, low_qual_spec, by = "sample", all.y = T)

high_qual_spec = fread("./annotations_vgp.txt")
high_qual_spec$sample = high_qual_spec[,c("Assc.")]
high_qual_spec_params = merge(fread("./estimated-spectra_4.txt"), fread("./estimated-parameters_4.txt"))
high_qual_spec_params$sample = gsub(pattern = ".hist", replacement = "", high_qual_spec_params$sample)
high_qual_spec_params_new = merge(high_qual_spec_params, high_qual_spec, by = "sample", all.y = T)
extend_lin = unique.data.frame(fread("./VGP_extended_lin.tsv")[,c("Family Scientific Name", "Extended lineage")])
high_qual_spec_params_new = merge(high_qual_spec_params_new, extend_lin, by.x = "Family", by.y = "Family Scientific Name", all.x = T, all.y = F, )

#write.csv(file = "./low_qual_genomes/curated_low_qual_species.tsv", x = low_qual_spec, quote = F, row.names = F)

low_qual_spec[!low_qual_spec$scientific_name %in% high_qual_spec$ScientificName, scientific_name]
tree = read.tree("./low_qual_genomes/roadies_v1.1.4_existing_genomes.tre") 
low_qual_spec[!gsub(" ", "_", low_qual_spec$scientific_name) %in% tree$tip.label, scientific_name]

high_qual_spec_params_new[,3:52] = (high_qual_spec_params_new[,3:52]/10^3)
high_qual_spec_params_new$genome_length = high_qual_spec_params_new$genome_length/1000

x = melt(data = high_qual_spec_params_new[,c("ScientificName", 
                                          "Superorder", "Extended lineage", "Lineage", "genome_length",
                                          "r1","r2","r3","r4","r5","r6","r7","r8","r9",
                                          "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
                                          "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
                                          "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
                                          "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50")], 
         id.vars = c("ScientificName", "Superorder", "Extended lineage", "Lineage", "genome_length"),
         measure.vars = c("r1","r2","r3","r4","r5","r6","r7","r8","r9",
         "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
         "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
         "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
         "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50"))
#x$spectra = "VGP"

low_qual_spec_params_new[,2:51] = (low_qual_spec_params_new[,2:51]/10^3)
#high_qual_spec_params_new$genome_length = high_qual_spec_params_new$genome_length/1000
y = melt(data = low_qual_spec_params_new[,c("scientific_name", 
                                             "r1","r2","r3","r4","r5","r6","r7","r8","r9",
                                             "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
                                             "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
                                             "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
                                             "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50")], 
         id.vars = c("scientific_name"),
         measure.vars = c("r1","r2","r3","r4","r5","r6","r7","r8","r9",
                          "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
                          "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
                          "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
                          "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50"))
names(y) = c("ScientificName", "variable", "value")
#y


z = merge(x=x, y=y, by = c("ScientificName", "variable"), all.y = T, all.x = F)[ScientificName != "Rhinophrynus dorsalis"]

z$Lineage  = factor(z$Lineage, levels = c("Birds",
                                          "Amphibians",
                                          "Fishes", 
                                          "Reptiles",
                                          "Mammals",
                                          "Ancestral"))


sign_y <- function(x){sign(x)*log(abs(x))}

lineage = ggplot(data = z[!is.na(`Extended lineage`)], 
                 aes(
                     x = as.numeric(sub("r", "", variable)),
                     #y = sin_y((value.x - value.y)/value.y),
                     y = (value.x - value.y)/value.y,
                     color = `Extended lineage`)) +
  facet_wrap(.~`Lineage`, nrow = 1) +
  #stat_smooth(method = "loess") +
  stat_summary(fun = median, geom = "line", size = 1.5, alpha = 0.7) +
  #stat_summary(fun.data = mean_se, geom = "errorbar") +
  #stat_summary(fun = median) +
  #geom_point(alpha = 0.5) +
  theme_classic() + 
  ylab("percent increase") + xlab("k-mer frequency bin") + 
  #scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = c("Mammals" = "#http://127.0.0.1:37257/graphics/plot_zoom_png?width=2745&height=563E69F00",  # Okabe–Ito Orange
  #                               "Birds" = "#00796B",    # Teal 700
  #                               "Crocodiles" = "#009688",# Teal 500
  #                               "Turtles" = "#4DB6AC",  # Teal 300
  #                               "Lepidosauria" = "#80CBC4", # Teal 200
  #                               "Amphibians" = "#984EA3",   # Purple
  #                               "Lobe-finned fishes" = "#A6761D", # Brown
  #                               "Ray-finned fishes" = "#56B4E9",  # Okabe–Ito Sky Blue
  #                               "Cartilaginous fishes" = "#0072B2", # Okabe–Ito Blue
  #                               "Cyclostomes" = "#CC79A7", # Okabe–Ito Magenta
  #                               "Ancestral" = "#999999"  # Neutral Gray
  # )) +   
  theme(axis.text.x = element_text(size = 15),
               axis.text.y = element_text(size=15),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
               title = element_text(size=15),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
               panel.grid.minor.x = element_blank(),
               panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
               legend.position = "bottom",
               legend.text = element_text(size=15),
               legend.title = element_text(size=20),
               legend.key.spacing.x = unit(2.0, "cm"),
               strip.text = element_text(size = 17),
               panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
               plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(transform = scales::transform_pseudo_log(),labels = scales::percent) 
  #scale_y_continuous(transform = scales::transform_pseudo_log()) + scale_x_continuous(transform = scales::transform_pseudo_log())
lineage
lineage = NULL

z$i=as.numeric(sub("r", "", z$variable))

head(z)
  library(dplyr)
library(scales)
arrow_plot = z[!is.na(`Extended lineage`)] %>%
  group_by(`Extended lineage`,i) %>%
  summarize(newm = mean(value.y*i/genome_length),oldm = mean(value.x*i/genome_length)) %>%
  filter(i>1) %>%
ggplot(
                 aes(
                   x = i,
                   #y = sin_y((value.x - value.y)/value.y),
                   yend =  newm,
                   y = oldm ,
                   color = `Extended lineage`)) +
  geom_segment(size = 0.65, alpha = 1,position = position_dodge(width=0.5),
               arrow = arrow(ends = "first",length =  unit(4,"pt"))) +
  facet_wrap(.~`Extended lineage`, nrow = 2) +
  
  #stat_summary(fun.data = mean_se, geom = "errorbar") +
  #stat_summary(fun = median) +
  #geom_point(alpha = 0.5) +
  theme_classic() + 
  scale_color_manual(name="", values = colorRampPalette(RColorBrewer::brewer.pal(n=8, name = "Accent"))(10)) +
  #scale_color_brewer(palette = "Dark2") +
   scale_color_manual(name="",
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
        legend.position = "bottom",
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(name="portion of genome",labels=percent,trans="log10")+
  scale_x_continuous(breaks=c(2,10,20,30,40,50),name="31-mer frequency")
arrow_plot
ggsave("./arrow_plot.pdf", plot = arrow_plot, width = 16, height = 7)
