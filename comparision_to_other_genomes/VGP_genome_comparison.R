library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)
library(dplyr)

setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/comparision_to_other_genomes/")
source("./subset_phylogeny.R")

# =========================================== #
# MAXIMUM LIKELIHOOD ANCESTRAL RECONSTRUCTION #
# =========================================== #
# ------------------- HIGH COPY REPEAT PER MILLION ----------------------- #
# ~ VGP ~
vgp_HCRM_annots = setNames(vgp_merged$HCRM, vgp_merged$ScientificName)
vgp_HCRM_ml_asr = anc.ML(tree, vgp_HCRM_annots, model = "BM")
writeAncestors(tree = tree, Anc = vgp_HCRM_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_VGP_HCRM-ASR.nwk")
# extract internal node states
vgp_HCRM_states = as.data.table(vgp_HCRM_ml_asr$ace, keep.rownames = T)
names(vgp_HCRM_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(vgp_merged, keep.rownames = T)[,c("ScientificName", "HCRM")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with VGP annotations
vgp_HCRM_states = rbind(vgp_HCRM_states, tip_labels, use.names = TRUE, fill = TRUE)
vgp_HCRM_states = merge(vgp_HCRM_states, as.data.table(vgp_merged, keep.rownames = T), all = T)
vgp_HCRM_states$node = as.numeric(vgp_HCRM_states$node)

# ~ Existing Genomes ~
existing_HCRM_annots = setNames(existing_merged$HCRM, existing_merged$ScientificName)
existing_HCRM_ml_asr = anc.ML(tree, existing_HCRM_annots, model = "BM")
writeAncestors(tree = tree, Anc = existing_HCRM_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_existing_HCRM-ASR.nwk")
# extract internal node states
existing_HCRM_states = as.data.table(existing_HCRM_ml_asr$ace, keep.rownames = T)
names(existing_HCRM_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(existing_merged, keep.rownames = T)[,c("ScientificName", "HCRM")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with Existing Genome annotations
existing_HCRM_states= rbind(existing_HCRM_states, tip_labels, use.names = TRUE, fill = TRUE)
existing_HCRM_states = merge(existing_HCRM_states, as.data.table(existing_merged, keep.rownames = T), all = T)
existing_HCRM_states$node = as.numeric(existing_HCRM_states$node)

# ------------------- UNIQUENESS RATIO ----------------------- #
# ~ VGP ~
vgp_UR_annots = setNames(vgp_merged$repeat_ratio, vgp_merged$ScientificName)
vgp_UR_ml_asr = anc.ML(tree, vgp_UR_annots, model = "BM")
writeAncestors(tree = tree, Anc = vgp_UR_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_VGP_UR-ASR.nwk")
# extract internal node states
vgp_UR_states = as.data.table(vgp_UR_ml_asr$ace, keep.rownames = T)
names(vgp_UR_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(vgp_merged, keep.rownames = T)[,c("ScientificName", "repeat_ratio")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with VGP annotations
vgp_UR_states = rbind(vgp_UR_states, tip_labels, use.names = TRUE, fill = TRUE)
vgp_UR_states = merge(vgp_UR_states, as.data.table(vgp_merged, keep.rownames = T), all = T)
vgp_UR_states$node = as.numeric(vgp_UR_states$node)

# ~ Existing Genomes ~
existing_UR_annots = setNames(existing_merged$repeat_ratio, existing_merged$ScientificName)
existing_UR_ml_asr = anc.ML(tree, existing_UR_annots, model = "BM")
writeAncestors(tree = tree, Anc = existing_UR_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_existing_UR-ASR.nwk")
# extract internal node states
existing_UR_states = as.data.table(existing_UR_ml_asr$ace, keep.rownames = T)
names(existing_UR_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(existing_merged, keep.rownames = T)[,c("ScientificName", "repeat_ratio")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with Existing Genome annotations
existing_UR_states= rbind(existing_UR_states, tip_labels, use.names = TRUE, fill = TRUE)
existing_UR_states = merge(existing_UR_states, as.data.table(existing_merged, keep.rownames = T), all = T)
existing_UR_states$node = as.numeric(existing_UR_states$node)

# ------------------- GENOME LENGTH ----------------------- #
# ~ VGP ~
vgp_GL_annots = setNames(vgp_merged$genome_length, vgp_merged$ScientificName)
vgp_GL_ml_asr = anc.ML(tree, vgp_GL_annots, model = "BM")
writeAncestors(tree = tree, Anc = vgp_GL_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_VGP_GL-ASR.nwk")
# extract internal node states
vgp_GL_states = as.data.table(vgp_GL_ml_asr$ace, keep.rownames = T)
names(vgp_GL_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(vgp_merged, keep.rownames = T)[,c("ScientificName", "genome_length")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with VGP annotations
vgp_GL_states = rbind(vgp_GL_states, tip_labels, use.names = TRUE, fill = TRUE)
vgp_GL_states = merge(vgp_GL_states, as.data.table(vgp_merged, keep.rownames = T), all = T)
vgp_GL_states$node = as.numeric(vgp_GL_states$node)

# ~ Existing Genomes ~
existing_GL_annots = setNames(existing_merged$genome_length, existing_merged$ScientificName)
existing_GL_ml_asr = anc.ML(tree, existing_GL_annots, model = "BM")
writeAncestors(tree = tree, Anc = existing_GL_ml_asr, file = "../comparision_to_other_genomes/ASR/subset_existing_GL-ASR.nwk")
# extract internal node states
existing_GL_states = as.data.table(existing_GL_ml_asr$ace, keep.rownames = T)
names(existing_GL_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(existing_merged, keep.rownames = T)[,c("ScientificName", "genome_length")], 
                   by.x = "V2", by.y = "ScientificName")
names(tip_labels) = c('ScientificName', 'node', 'state')
# merging ASR data with Existing Genome annotations
existing_GL_states= rbind(existing_GL_states, tip_labels, use.names = TRUE, fill = TRUE)
existing_GL_states = merge(existing_GL_states, as.data.table(existing_merged, keep.rownames = T), all = T)
existing_GL_states$node = as.numeric(existing_GL_states$node)


# =========================================== #
#        RECONSTRUCTION SCATTER PLOT          #
# =========================================== #

# ------------------- HCRM ----------------------- #
all_HCRM_states = merge(x = vgp_HCRM_states, y = existing_HCRM_states, by = "node", all = T)
all_HCRM_states[is.na(all_HCRM_states$ScientificName.x), `Extended lineage` := "Ancestral"]
all_HCRM_states$ancestral = ifelse(all_HCRM_states$`Extended lineage` == "Ancestral", T, F)
all_HCRM_states$`Extended lineage` = factor(all_HCRM_states$`Extended lineage`, levels = c("Amphibians", "Birds", "Cartilaginous fishes",
                                                                                       "Crocodiles", "Cyclostomes", "Lepidosauria", 
                                                                                       "Lobe-finned fishes", "Mammals", "Ray-finned fishes",
                                                                                       "Turtles", "Ancestral"))

#SCATTER PLOT
HCRM_plot = ggplot(data = all_HCRM_states[!is.na(ScientificName.x)], 
                   aes(x = state.y, y = state.x, color = `Extended lineage`, shape = ancestral)) + 
  geom_abline() +
  geom_point(data = all_HCRM_states[is.na(ScientificName.x)], alpha=0.7) +
  geom_point(alpha=0.7) +
  theme_classic() +
  guides(shape = "none", color = guide_legend(reverse = T)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        title = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
        legend.position = "left",
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_shape_manual(values = c(19,17)) +
  xlab("Existing Genomes") + 
  ylab("") + 
  ggtitle("HCRM") +
  guides(color = "none") +  
  #guides(color = guide_legend("Lineage")) + 
  scale_color_manual(values = c("Mammals" = "#E69F00",  # Okabe–Ito Orange
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
  coord_cartesian(xlim=c(0,700), ylim=c(0,700))
#HCRM_plot


# ------------------- UNIQUENESS RATIO ----------------------- #
all_UR_states = merge(x = vgp_UR_states, y = existing_UR_states, by = "node", all = T)
all_UR_states[is.na(all_UR_states$ScientificName.x), `Extended lineage` := "Ancestral"]
all_UR_states$ancestral = ifelse(all_UR_states$`Extended lineage` == "Ancestral", T, F)
all_UR_states$`Extended lineage` = factor(all_UR_states$`Extended lineage`, levels = c("Amphibians", "Birds", "Cartilaginous fishes",
                                                                                           "Crocodiles", "Cyclostomes", "Lepidosauria", 
                                                                                           "Lobe-finned fishes", "Mammals", "Ray-finned fishes",
                                                                                           "Turtles", "Ancestral"))

#SCATTER PLOT
UR_plot = ggplot(data = all_UR_states[!is.na(ScientificName.x)], 
                 aes(x = state.y, state.x, color = `Extended lineage`, shape = `Extended lineage`)) + 
  geom_abline() +
  geom_point(data = all_UR_states[is.na(ScientificName.x)], alpha=0.7, size = 1.7) +
  geom_point(alpha=0.7, size = 1.7) +
  theme_classic() +
  guides(shape = "none", color = guide_legend(reverse = T)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        title = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
        legend.position = "left",
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_shape_manual(values = c(19,17)) +
  xlab("") + 
  ylab("VGP Genomes") + 
  ggtitle("Repetitive 31-mers") +
  guides(color = guide_legend("Lineage", override.aes = list(size = 3)), shape = guide_legend("Lineage")) + 
  scale_shape_manual(values = c("Mammals" = 16,  # Okabe–Ito Orange
                                "Birds" = 16,    # Teal 700
                                "Crocodiles" = 16,# Teal 500
                                "Turtles" = 16,  # Teal 300
                                "Lepidosauria" = 16, # Teal 200
                                "Amphibians" = 16,   # Purple
                                "Lobe-finned fishes" = 16, # Brown
                                "Ray-finned fishes" = 16,  # Okabe–Ito Sky Blue
                                "Cartilaginous fishes" = 16, # Okabe–Ito Blue
                                "Cyclostomes" = 16, # Okabe–Ito Magenta
                                "Ancestral" = 17  # Neutral Gray
  )) +
  scale_color_manual(values = c("Mammals" = "#E69F00",  # Okabe–Ito Orange
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
  #scale_color_manual(values = c("grey50", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) +
  coord_cartesian(xlim=c(0,0.8), ylim=c(0,0.8))
#UR_plot

#"#A50026" "#D73027" "#F46D43" "#FDAE61" "#FEE090" "#E0F3F8" "#ABD9E9" "#74ADD1" "#4575B4" "#313695"
#COMPARISON TREE

# ------------------- GENOME LENGTH ----------------------- #
all_GL_states = merge(x = vgp_GL_states, y = existing_GL_states, by = "node", all = T)
all_GL_states$state.x = as.numeric(all_GL_states$state.x)/1000000
all_GL_states$state.y = as.numeric(all_GL_states$state.y)/1000000
all_GL_states[is.na(all_GL_states$ScientificName.x), `Extended lineage` := "Ancestral"]
all_GL_states$ancestral = ifelse(all_GL_states$`Extended lineage` == "Ancestral", T, F)
all_GL_states$`Extended lineage` = factor(all_GL_states$`Extended lineage`, levels = c("Amphibians", "Birds", "Cartilaginous fishes",
                                                                                           "Crocodiles", "Cyclostomes", "Lepidosauria", 
                                                                                           "Lobe-finned fishes", "Mammals", "Ray-finned fishes",
                                                                                           "Turtles", "Ancestral"))
#SCATTER PLOT
GL_plot = ggplot(data = all_GL_states[!is.na(ScientificName.x)], 
       aes(x = state.y, state.x, color = `Extended lineage`, shape = ancestral)) + 
  geom_abline() +
  geom_point(data = all_GL_states[is.na(ScientificName.x)], alpha=0.7) +
  geom_point(alpha=0.7) +
  theme_classic() +
  guides(shape = "none", color = guide_legend(reverse = T)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        title = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
        legend.position = "left",
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_shape_manual(values = c(19,17)) +
  xlab("") + 
  ylab("") + 
  ggtitle("Genome Length (Mbp)") +
  guides(color = "none") + 
  scale_color_manual(values = c("Mammals" = "#E69F00",  # Okabe–Ito Orange
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
  coord_cartesian(xlim=c(0,4000), ylim=c(0,4000))
#GL_plot

ggsave(plot = UR_plot + HCRM_plot + GL_plot, "../figures/all_plots.pdf", width = 25*0.8 , height = 7*0.8)


####################
# COMPARISON TREES #
####################

scale_percent = function(x){scales::percent(x, scale=1)}

####################
#HCRM_comp_tree

HCRM_comp_tree = ggtree(tree, color = "black", size = 1.5, branch.length="none") %<+% all_HCRM_states +
  layout_circular() +  
  geom_tree(aes(color = (state.x-state.y)/state.y * 100 ),
            continuous = 'colour',
            size=1.3) +
  scale_color_gradient2(low = "#d500d5",
                        mid = "#ffffff",
                        high = "#008020",
                        midpoint = 0,
                        transform = "pseudo_log",
                        breaks = c(-10, 0, 10, 100, 1000, 10000),
                        labels = scale_percent) +
  guides(color = guide_colorbar("% Change in HCRM")) +
  geom_strip(taxa1 = "Myxine_glutinosa", taxa2 = "Lampetra_fluviatilis", extend = 0.3,
             label = "Agnatha", offset.text=1, angle = 1, 
             color = "#000") + 
  geom_strip(taxa1 = "Stegostoma_tigrinum", taxa2 = "Carcharodon_carcharias",extend = 0.3,
             label = "Chondrichthyes", offset.text=1, angle = 5, 
             color = "#000") +
  geom_strip(taxa1 = "Acipenser_ruthenus", taxa2 = "Sparus_aurata", extend = 0.3,
             label = "Actinopterygii", offset.text=1, angle = 55,
             color = "#000") + 
  geom_strip(taxa1 = "Protopterus_annectens", taxa2 = "Latimeria_chalumnae", extend = 0.3,
             label = "Sarcopterygii", offset.text=8, angle = 280,
             color = "#000") +
  geom_strip(taxa1 = "Ascaphus_truei", taxa2 = "Engystomops_pustulosus",extend = 0.3,
             label = "Anura", offset.text=4.5, angle = 288,
             color = "#000") +
  geom_strip(taxa1 = "Ornithorhynchus_anatinus", taxa2 = "Tachyglossus_aculeatus",extend = 0.3,
             label = "Monotremes", offset.text=7.5, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Monodelphis_domestica", taxa2 = "Sarcophilus_harrisii",extend = 0.3,
             label = "Marsupials", offset.text=7, angle = 300,
             color = "#000") +
  geom_strip(taxa1 = "Dasypus_novemcinctus", taxa2 = "Tamandua_tetradactyla",extend = 0.3,
             label = "Xenarthra", offset.text=6.5, angle = 304,
             color = "#000") +
  geom_strip(taxa1 = "Dugong_dugon", taxa2 = "Loxodonta_africana",extend = 0.3,
             label = "Afrotheria", offset.text=6.5, angle = 310,
             color = "#000") +
  geom_strip(taxa1 = "Lemur_catta", taxa2 = "Urocitellus_parryii",extend = 0.3,
             label = "Supraprimates", offset.text=9, angle = 330,
             color = "#000") +
  geom_strip(taxa1 = "Sorex_araneus", taxa2 = "Grampus_griseus",extend = 0.3,
             label = "Laurasiatheria", offset.text=9, angle = 25,
             color = "#000") +
  geom_strip(taxa1 = "Tiliqua_scincoides", taxa2 = "Podarcis_siculus",extend = 0.3,
             label = "Squamata", offset.text=6.5, angle = 60, 
             color = "#000") +
  geom_strip(taxa1 = "Podocnemis_expansa", taxa2 = "Podocnemis_expansa",
             label = "\nTestudines (Pleurodira)", offset.text=13, angle = 65, extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Carettochelys_insculpta", taxa2 = "Chelonia_mydas",
             label = "\nTestudines (Cryptodira)", offset.text=13, angle = 65,extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Gavialis_gangeticus", taxa2 = "Alligator_mississippiensis",extend = 0.3,
             label = "Crocodylia", offset.text=6.5, angle = 70,
             color = "#000") +
  geom_strip(taxa1 = "Eudromia_elegans", taxa2 = "Rhea_pennata",extend = 0.3,
             label = "Palaeognathae", offset.text = 8.5, angle = 260 + 170,
             color = "#000") +
  geom_strip(taxa1 = "Gallus_gallus", taxa2 = "Anser_cygnoides",extend = 0.3,
             label = "GalloAnserformes", offset.text = 10, angle = 260+180,
             color = "#000") +
  geom_strip(taxa1 = "Cuculus_canorus", taxa2 = "Chlamydotis_macqueenii",extend = 0.3,
             label = "Otidimorphae", offset.text = 1, angle = 273,
             color = "#000") +
  geom_strip(taxa1 = "Phoenicopterus_ruber", taxa2 = "Caloenas_nicobarica",extend = 0.3,
             label = "Columbea", offset.text = 1, angle = 275,
             color = "#000") +
  geom_strip(taxa1 = "Caprimulgus_europaeus", taxa2 = "Podargus_strigoides",extend = 0.3,
             label = "Caprimulgimorphae", offset.text = 1, angle = 270,
             color = "#000") +
  geom_strip(taxa1 = "Opisthocomus_hoazin", taxa2 = "Larus_michahellis",extend = 0.3,
             label = "Core waterbirds", offset.text = 1, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Strix_aluco", taxa2 = "Agelaius_phoeniceus",extend = 0.3,
             label = "Core landbirds", offset.text = 1, angle = 332,
             color = "#000") +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm")) 
#ggsave(filename = "./HCRM_comp_tree.pdf", plot = HCRM_comp_tree, height = 14, width = 15)

## UR COMP TREE
#summary(with(all_UR_states,(state.x-state.y)/state.y ))

UR_comp_tree = ggtree(tree, color = "black", size = 1.5, branch.length="none") %<+% all_UR_states +
  layout_circular() +  
  geom_tree(aes(color = (state.x-state.y)/state.y * 10 ), 
            #continuous = 'colour', 
            size=1.3) +
  guides(color = guide_colorbar("% Change in RR")) +
  geom_strip(taxa1 = "Myxine_glutinosa", taxa2 = "Lampetra_fluviatilis", extend = 0.3,
             label = "Agnatha", offset.text=1, angle = 1, 
             color = "#000") + 
  geom_strip(taxa1 = "Stegostoma_tigrinum", taxa2 = "Carcharodon_carcharias",extend = 0.3,
             label = "Chondrichthyes", offset.text=1, angle = 5, 
             color = "#000") +
  geom_strip(taxa1 = "Acipenser_ruthenus", taxa2 = "Sparus_aurata", extend = 0.3,
             label = "Actinopterygii", offset.text=1, angle = 55,
             color = "#000") + 
  geom_strip(taxa1 = "Protopterus_annectens", taxa2 = "Latimeria_chalumnae", extend = 0.3,
             label = "Sarcopterygii", offset.text=8, angle = 280,
             color = "#000") +
  geom_strip(taxa1 = "Ascaphus_truei", taxa2 = "Engystomops_pustulosus",extend = 0.3,
             label = "Anura", offset.text=4.5, angle = 288,
             color = "#000") +
  geom_strip(taxa1 = "Ornithorhynchus_anatinus", taxa2 = "Tachyglossus_aculeatus",extend = 0.3,
             label = "Monotremes", offset.text=7.5, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Monodelphis_domestica", taxa2 = "Sarcophilus_harrisii",extend = 0.3,
             label = "Marsupials", offset.text=7, angle = 300,
             color = "#000") +
  geom_strip(taxa1 = "Dasypus_novemcinctus", taxa2 = "Tamandua_tetradactyla",extend = 0.3,
             label = "Xenarthra", offset.text=6.5, angle = 304,
             color = "#000") +
  geom_strip(taxa1 = "Dugong_dugon", taxa2 = "Loxodonta_africana",extend = 0.3,
             label = "Afrotheria", offset.text=6.5, angle = 310,
             color = "#000") +
  geom_strip(taxa1 = "Lemur_catta", taxa2 = "Urocitellus_parryii",extend = 0.3,
             label = "Supraprimates", offset.text=9, angle = 330,
             color = "#000") +
  geom_strip(taxa1 = "Sorex_araneus", taxa2 = "Grampus_griseus",extend = 0.3,
             label = "Laurasiatheria", offset.text=9, angle = 25,
             color = "#000") +
  geom_strip(taxa1 = "Tiliqua_scincoides", taxa2 = "Podarcis_siculus",extend = 0.3,
             label = "Squamata", offset.text=6.5, angle = 60, 
             color = "#000") +
  geom_strip(taxa1 = "Podocnemis_expansa", taxa2 = "Podocnemis_expansa",
             label = "\nTestudines (Pleurodira)", offset.text=13, angle = 65, extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Carettochelys_insculpta", taxa2 = "Chelonia_mydas",
             label = "\nTestudines (Cryptodira)", offset.text=13, angle = 65,extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Gavialis_gangeticus", taxa2 = "Alligator_mississippiensis",extend = 0.3,
             label = "Crocodylia", offset.text=6.5, angle = 70,
             color = "#000") +
  geom_strip(taxa1 = "Eudromia_elegans", taxa2 = "Rhea_pennata",extend = 0.3,
             label = "Palaeognathae", offset.text = 8.5, angle = 260 + 170,
             color = "#000") +
  geom_strip(taxa1 = "Gallus_gallus", taxa2 = "Anser_cygnoides",extend = 0.3,
             label = "GalloAnserformes", offset.text = 10, angle = 260+180,
             color = "#000") +
  geom_strip(taxa1 = "Cuculus_canorus", taxa2 = "Chlamydotis_macqueenii",extend = 0.3,
             label = "Otidimorphae", offset.text = 1, angle = 273,
             color = "#000") +
  geom_strip(taxa1 = "Phoenicopterus_ruber", taxa2 = "Caloenas_nicobarica",extend = 0.3,
             label = "Columbea", offset.text = 1, angle = 275,
             color = "#000") +
  geom_strip(taxa1 = "Caprimulgus_europaeus", taxa2 = "Podargus_strigoides",extend = 0.3,
             label = "Caprimulgimorphae", offset.text = 1, angle = 270,
             color = "#000") +
  geom_strip(taxa1 = "Opisthocomus_hoazin", taxa2 = "Larus_michahellis",extend = 0.3,
             label = "Core waterbirds", offset.text = 1, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Strix_aluco", taxa2 = "Agelaius_phoeniceus",extend = 0.3,
             label = "Core landbirds", offset.text = 1, angle = 332,
             color = "#000") +
  scale_color_gradient2(low = "#d500d5",
                        mid = "#ffffff",
                        high = "#008020",
                        midpoint = 0,
                        transform = "pseudo_log",
                        breaks = c(-30, 0, 30, 100, 350, 1000)/10,
                        labels = function(x){scales::percent(x, scale=10)}) +
  # scale_color_stepsn(colors = rev(c("#D73027", "#FDAE61","#FEF090","#FFFFFF", "#ABD9E9")),
  #                    nice.breaks = T,
  #                    transform = "pseudo_log",
  #                    breaks = c(-75,-25, 0, 25, 100, 300,1000)/10,
  #                    labels = function(x){scales::percent(x, scale=10)}) +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm")) 

#ggsave(filename = "../figures/RR_comparison_tree.pdf", plot = UR_comp_tree, height = 14, width = 15)



#COMPARISON TREE
GL_comp_tree = ggtree(tree, color = "black", size = 1.5, branch.length="none") %<+% all_GL_states +
  layout_circular() +   
  geom_tree(aes(color = (state.x-state.y)/state.y * 100 ), 
            continuous = 'colour', 
            size=1.3) +
  scale_color_gradient2(low = "#d500d5", 
                        mid = "#ffffff", 
                        high = "#008020", 
                        midpoint = 0, 
                        transform = "pseudo_log", 
                        breaks = c(-10, 0, 10, 100), 
                        labels = scale_percent) +
  guides(color = guide_colorbar("% Change in Genome Length")) +
  geom_strip(taxa1 = "Myxine_glutinosa", taxa2 = "Lampetra_fluviatilis", extend = 0.3,
             label = "Agnatha", offset.text=1, angle = 1, 
             color = "#000") + 
  geom_strip(taxa1 = "Stegostoma_tigrinum", taxa2 = "Carcharodon_carcharias",extend = 0.3,
             label = "Chondrichthyes", offset.text=1, angle = 5, 
             color = "#000") +
  geom_strip(taxa1 = "Acipenser_ruthenus", taxa2 = "Sparus_aurata", extend = 0.3,
             label = "Actinopterygii", offset.text=1, angle = 55,
             color = "#000") + 
  geom_strip(taxa1 = "Protopterus_annectens", taxa2 = "Latimeria_chalumnae", extend = 0.3,
             label = "Sarcopterygii", offset.text=8, angle = 280,
             color = "#000") +
  geom_strip(taxa1 = "Ascaphus_truei", taxa2 = "Engystomops_pustulosus",extend = 0.3,
             label = "Anura", offset.text=4.5, angle = 288,
             color = "#000") +
  geom_strip(taxa1 = "Ornithorhynchus_anatinus", taxa2 = "Tachyglossus_aculeatus",extend = 0.3,
             label = "Monotremes", offset.text=7.5, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Monodelphis_domestica", taxa2 = "Sarcophilus_harrisii",extend = 0.3,
             label = "Marsupials", offset.text=7, angle = 300,
             color = "#000") +
  geom_strip(taxa1 = "Dasypus_novemcinctus", taxa2 = "Tamandua_tetradactyla",extend = 0.3,
             label = "Xenarthra", offset.text=6.5, angle = 304,
             color = "#000") +
  geom_strip(taxa1 = "Dugong_dugon", taxa2 = "Loxodonta_africana",extend = 0.3,
             label = "Afrotheria", offset.text=6.5, angle = 310,
             color = "#000") +
  geom_strip(taxa1 = "Lemur_catta", taxa2 = "Urocitellus_parryii",extend = 0.3,
             label = "Supraprimates", offset.text=9, angle = 330,
             color = "#000") +
  geom_strip(taxa1 = "Sorex_araneus", taxa2 = "Grampus_griseus",extend = 0.3,
             label = "Laurasiatheria", offset.text=9, angle = 25,
             color = "#000") +
  geom_strip(taxa1 = "Tiliqua_scincoides", taxa2 = "Podarcis_siculus",extend = 0.3,
             label = "Squamata", offset.text=6.5, angle = 60, 
             color = "#000") +
  geom_strip(taxa1 = "Podocnemis_expansa", taxa2 = "Podocnemis_expansa",
             label = "\nTestudines (Pleurodira)", offset.text=13, angle = 65, extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Carettochelys_insculpta", taxa2 = "Chelonia_mydas",
             label = "\nTestudines (Cryptodira)", offset.text=13, angle = 65,extend = 0.3,
             color = "#000") +
  geom_strip(taxa1 = "Gavialis_gangeticus", taxa2 = "Alligator_mississippiensis",extend = 0.3,
             label = "Crocodylia", offset.text=6.5, angle = 70,
             color = "#000") +
  geom_strip(taxa1 = "Eudromia_elegans", taxa2 = "Rhea_pennata",extend = 0.3,
             label = "Palaeognathae", offset.text = 8.5, angle = 260 + 170,
             color = "#000") +
  geom_strip(taxa1 = "Gallus_gallus", taxa2 = "Anser_cygnoides",extend = 0.3,
             label = "GalloAnserformes", offset.text = 10, angle = 260+180,
             color = "#000") +
  geom_strip(taxa1 = "Cuculus_canorus", taxa2 = "Chlamydotis_macqueenii",extend = 0.3,
             label = "Otidimorphae", offset.text = 1, angle = 273,
             color = "#000") +
  geom_strip(taxa1 = "Phoenicopterus_ruber", taxa2 = "Caloenas_nicobarica",extend = 0.3,
             label = "Columbea", offset.text = 1, angle = 275,
             color = "#000") +
  geom_strip(taxa1 = "Caprimulgus_europaeus", taxa2 = "Podargus_strigoides",extend = 0.3,
             label = "Caprimulgimorphae", offset.text = 1, angle = 270,
             color = "#000") +
  geom_strip(taxa1 = "Opisthocomus_hoazin", taxa2 = "Larus_michahellis",extend = 0.3,
             label = "Core waterbirds", offset.text = 1, angle = 295,
             color = "#000") +
  geom_strip(taxa1 = "Strix_aluco", taxa2 = "Agelaius_phoeniceus",extend = 0.3,
             label = "Core landbirds", offset.text = 1, angle = 332,
             color = "#000") +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm")) 


#UR_comp_tree + HCRM_comp_tree + GL_comp_tree
ggsave(plot = UR_comp_tree + HCRM_comp_tree + GL_comp_tree, "../figures/ASR_comparison_trees.pdf", width = 45 , height = 15)


########################
#mean and average increases in genome length
all_GL_states %>% 
  filter(ScientificName.x %in% tree$tip.label) %>% 
  summarize(vgp_avg = mean((state.x-state.y)/state.y), 
            vgp_med = median((state.x-state.y)/state.y))

#mean and average increases in repeat ratio
all_UR_states %>% 
  filter(ScientificName.x %in% tree$tip.label) %>% 
  summarize(vgp_avg = mean((state.x-state.y)/state.y), 
            vgp_med = median((state.x-state.y)/state.y))

#number of significantly less repetitive samples
all_UR_states %>% 
  filter(ScientificName.x %in% tree$tip.label) %>% 
  filter((state.x - state.y) / state.y <= -0.05) %>%
  nrow()

#mean and average increases in HCRM
all_HCRM_states %>% 
  filter(ScientificName.x %in% tree$tip.label) %>% 
  summarize(vgp_avg = mean((state.x-state.y)/state.y), 
            vgp_med = median((state.x-state.y)/state.y))

#Chimpanzee Stats
all_GL_states %>% 
  filter(ScientificName.x == "Pan_troglodytes") %>%
  summarize((state.x-state.y)/state.y, state.x, state.y)

all_UR_states %>% 
  filter(ScientificName.x == "Pan_troglodytes") %>%
  summarize((state.x-state.y)/state.y, state.x, state.y)  

#Bonobo Stats
all_GL_states %>% 
  filter(ScientificName.x == "Pan_paniscus") %>%
  summarize((state.x-state.y)/state.y, state.x, state.y)

all_UR_states %>% 
  filter(ScientificName.x == "Pan_paniscus") %>%
  summarize((state.x-state.y)/state.y, state.x, state.y)  

# GL Change from Bonobo to Chimp
(3239.869 - 3172.755)/ 3172.755 # = 2%
# RR Change from Bonobo to Chimp
(0.23-0.25)/0.25 # 8%
