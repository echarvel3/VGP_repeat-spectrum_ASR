
library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)

setwd("~/Desktop/VGP_analyses/")

# functions 
split_row <- function(x){
  z = paste(collapse = "_", strsplit(x, split = "_")[[1]][1:2])
}

sign_y <- function(x){sign(x)*log(abs(x))}

# reading data
low_qual_spec = fread("./low_qual_genomes/vertebrate_querySpeciesRep.tsv")
low_qual_spec$sample = low_qual_spec$assembly_id
low_qual_spec = low_qual_spec[,!c("assembly_id")]

# relabeling species
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

# reading spectrum data for existing genomes
low_qual_spec_params = merge(
  fread("./VGP_repeat-spectrum_ASR/comparision_to_other_genomes/respect_full-existing_spectra.txt"),
  fread("./VGP_repeat-spectrum_ASR/comparision_to_other_genomes/respect_full-existing_parameters.txt")
)
low_qual_spec_params$sample = sapply(X = low_qual_spec_params$sample, FUN = split_row)
low_qual_spec_params_new = merge(low_qual_spec_params, low_qual_spec, by = "sample", all.y = T)

# reading spectrum data for VGP genomes
high_qual_spec = fread("./annotations_vgp.txt")
high_qual_spec$sample = high_qual_spec[,c("Assc.")]

high_qual_spec_params = merge(
  fread("./estimated-spectra_4.txt"), 
  fread("./estimated-parameters_4.txt")
  )
high_qual_spec_params$sample = gsub(pattern = ".hist", replacement = "", high_qual_spec_params$sample)
high_qual_spec_params_new = merge(high_qual_spec_params, high_qual_spec, by = "sample", all.y = T)
extend_lin = unique.data.frame(fread("./VGP_extended_lin.tsv")[,c("Family Scientific Name", "Extended lineage")])
high_qual_spec_params_new = merge(high_qual_spec_params_new, extend_lin, by.x = "Family", by.y = "Family Scientific Name", all.x = T, all.y = F, )

# trimming data
tree = read.tree("./low_qual_genomes/roadies_v1.1.4_existing_genomes.tre") 
tree = drop.tip(tree, tip = c("Sminthopsis_crassicaudata",
                              "Shinisaurus_crocodilurus",
                              "Scardinius_erythrophthalmus",
                              "Tupaia_tana",
                              "Macrotis_lagotis",
                              "Aulostomus_maculatus",
                              "Salvelinus_alpinus",
                              "Spinachia_spinachia",
                              "Oenanthe_melanoleuca")
                )
low_qual_spec[!gsub(" ", "_", low_qual_spec$scientific_name) %in% tree$tip.label, scientific_name]


# formatting vgp data
high_qual_spec_params_new[,3:52] = (high_qual_spec_params_new[,3:52]/1000)
high_qual_spec_params_new$genome_length = high_qual_spec_params_new$genome_length/1000

x = melt(data = high_qual_spec_params_new[,c("ScientificName", 
                                          "Superorder", "Extended lineage", "Lineage", "genome_length",
                                          "r1","r2","r3","r4","r5","r6","r7","r8","r9",
                                          "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
                                          "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
                                          "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
                                          "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50")], 
         id.vars = c("ScientificName", "Superorder", "Extended lineage", "Lineage", "genome_length"),
         measure.vars = c(6:55))


# formatting existing data
low_qual_spec_params_new[,2:51] = (low_qual_spec_params_new[,2:51]/1000)
low_qual_spec_params_new$genome_length = as.numeric(low_qual_spec_params_new$genome_length)/1000
y = melt(data = low_qual_spec_params_new[,c("scientific_name", "genome_length", 
                                             "r1","r2","r3","r4","r5","r6","r7","r8","r9",
                                             "r10","r11","r12","r13","r14","r15","r16","r17","r18","r19",
                                             "r20","r21","r22","r23","r24","r25","r26","r27","r28","r29",
                                             "r30","r31","r32","r33","r34","r35","r36","r37","r38","r39",
                                             "r40","r41","r42","r43","r44","r45","r46","r47","r48","r49", "r50")], 
         id.vars = c("scientific_name", "genome_length"),
         measure.vars = c(3:52))
names(y) = c("ScientificName", "genome_length", "variable", "value")

z = merge(x=x, y=y, by = c("ScientificName", "variable"), all.y = T, all.x = F)[ScientificName != "Rhinophrynus dorsalis"]
z$i = as.numeric(gsub("r", "", z$variable))

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
        legend.position = c(0.9215, 0.79),
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
ggsave("./arrow_plot.pdf", plot = arrow_plot, width = 16, height = 7)


 z[!is.na(`Extended lineage`)] %>%
  arrange(`ScientificName`,-i) %>%
  group_by(`ScientificName`) %>% #print(n=70) 
  mutate(newm = as.numeric(cumsum(value.y*i)),
            oldm = as.numeric(cumsum(value.x*i))) %>% #print(n=70) 
   group_by(`Extended lineage`,i) %>%
   summarize(newmcc = mean(newm),oldmcc = mean(oldm)) %>% #print(n=70) 
ggplot(
  aes(
    x = i,
    #y = sin_y((value.x - value.y)/value.y),
    yend =  as.numeric(newmcc*1000),
    y = as.numeric(oldmcc*1000),
    color = `Extended lineage`)) +
  geom_segment(size = 0.65, alpha = 1,position = position_dodge(width=0.5),
               arrow = arrow(ends = "first",length =  unit(4,"pt"))) +
   facet_wrap(.~ifelse(
     `Extended lineage` %in% c("Ray-finned fishes", "Cyclostomes", "Birds","Lepidosauria"), 
     "Small", 
     ifelse(`Extended lineage` %in% c("Lobe-finned fishes"), "Large", "Medium")), 
     # "Ray-finned fishes", "Cartilaginous fishes", "Cyclostomes"
     nrow = 1, scales = "free_y") +
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
        legend.position =  "bottom",
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
