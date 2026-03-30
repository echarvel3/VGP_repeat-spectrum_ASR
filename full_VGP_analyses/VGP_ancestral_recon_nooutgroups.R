library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)


setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/full_VGP_analyses/")
#Read Full VGP Tree
tree = read.tree("../roadies_v1.1.16b.nwk") 
#Read Taxonomy Annotations
merged_annots = fread("../annotations_vgp.txt")
unique(merged_annots$VGPLineage)
tree = drop.tip(tree, merged_annots[Lineage %in% c("Invertebrates"), Assc.])
merged_annots = merged_annots[!VGPLineage %in% c("Invertebrates") & Assc. %in% tree$tip.label]
#Read VGP Respect Results
respect_results = fread("./respect_full-VGP_parameters.txt")
respect_results$sample = sub(".hist", "", respect_results$sample)
#Pool Data
merged_annots  = merge(x = respect_results,  y = merged_annots, by.x = "sample", by.y = "Assc.",)

# =========================================== #
# MAXIMUM LIKELIHOOD ANCESTRAL RECONSTRUCTION #
# =========================================== #

# Uniqueness Ratio ASR
UR_annots = setNames(object = 1-merged_annots$uniqueness_ratio, nm = merged_annots$sample)
UR_ml_asr = anc.ML(tree, UR_annots, model = "BM")
# writeAncestors(tree = tree, Anc = UR_ml_asr, file = "./full_VGP_UR-ASR.nwk")
# extract internal node states
UR_states = as.data.table(UR_ml_asr$ace, keep.rownames = T)
names(UR_states) = c("node", "state")
# adding leaf node states
tip_labels = as.data.table(cbind(c(1:length(tree$tip.label)), c(tree$tip.label)))
tip_labels = merge(tip_labels, as.data.table(UR_annots, keep.rownames = T), by.x = "V2", by.y = "rn")
names(tip_labels) = c('sample', 'node', 'state')
# merging ASR data with VGP annotations
UR_states = rbind(UR_states, tip_labels, use.names = T, fill = T)
UR_states = merge(UR_states, as.data.table(merged_annots, keep.rownames = T), all = T)
UR_states$node = as.numeric(UR_states$node)


# =========================================== #
#             TREE VISUALIZATION              #
# =========================================== #

# Uniqueness Ratio Tree
UR_tree_bl = ggtree(tree, 
                 color = "black", 
                 size = 0, 
                 #branch.length="none"
                 ) %<+% UR_states +
  geom_tree(aes(color = state), continuous = 'colour', size=1) +
  scale_color_gradientn(colours=rev(c("red", 'orange', 'green', 'cyan', 'blue')),
                        breaks=c(0.2, 0.4, 0.6),
                        limits=c(0, 0.8))+
  layout_circular() + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(24)) +
  guides(color = "none") +
  # geom_strip(taxa1 = "GCF_902459465.1", taxa2 = "GCF_902459465.1",
  #            label = "Echinodermata", offset.text = 0, angle = 0,
  #            color = "#000", extend = 0.3) +
  # geom_strip(taxa1 = "GCA_040954625.2", taxa2 = "GCA_040954625.2",
  #            label = "Hemichordata", offset.text = 0, angle = 2,
  #            color = "#000", extend = 0.3) +
  # geom_strip(taxa1 = "GCA_015852565.1", taxa2 = "GCA_019207075.1",
  #            label = "Cephalochordata", offset.text = , angle = 3,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_918807975.1", taxa2 = "GCA_965202585.1",
  #            label = "Tunicata", offset.text = 1, angle = 5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_964187855.1", taxa2 = "GCA_964198595.1",
  #            label = "Agnatha", offset.text=1, angle = 1+4,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_035084275.1", taxa2 = "GCA_964213995.1",
  #            label = "Chondrichthyes", offset.text=1, angle = 1+8,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_900747795.2", taxa2 = "GCF_904848185.1",
  #            label = "Actinopterygii", offset.text=1, angle = 60+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_040939525.1", taxa2 = "GCF_037176945.1",
  #            label = "Sarcopterygii", offset.text=10, angle = 289.5+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_901001135.1", taxa2 = "GCF_902459505.1",
  #            label = "Gymnophiona", offset.text=11, angle = 291+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_040938575.1", taxa2 = "GCA_964261635.1",
  #            label = "Caudata", offset.text=7, angle = 291.5+7,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_040206685.1", taxa2 = "GCF_040894005.1",
  #            label = "Anura", offset.text=5, angle = 300+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_004115215.2", taxa2 = "GCF_015852505.1",
  #            label = "Monotremes", offset.text=9, angle = 310+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_027887165.1", taxa2 = "GCF_011100635.1",
  #            label = "Marsupials", offset.text=8, angle = 315+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_030445035.2", taxa2 = "GCA_023851605.1",
  #            label = "Xenartha", offset.text=7, angle = 319+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_043290085.1", taxa2 = "GCF_030014295.1",
  #            label = "Afrotheria", offset.text=7, angle = 320+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_026018925.1", taxa2 = "GCA_949316315.1",
  #            label = "Supraprimates", offset.text=11, angle = 340+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_964194135.1", taxa2 = "GCF_011762595.1",
  #            label = "Laurasiatheria", offset.text=11, angle = 30+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_046126795.1", taxa2 = "GCA_964234715.1",
  #            label = "Squamata", offset.text=8, angle = 65+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_045364815.1", taxa2 = "GCA_040894355.2",
  #            label = "Testudines (Pleurodira)", offset.text=16, angle = 79+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_033958435.1", taxa2 = "GCA_965140285.1",
  #            label = "Testudines (Cryptodira)", offset.text=16, angle = 87+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_030020295.1", taxa2 = "GCF_030867095.1",
  #            label = "Crocodylia", offset.text = 1, angle = 270+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_040807025.1", taxa2 = "GCF_036370855.1",
  #            label = "Palaeognathae", offset.text = 1, angle = 275+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_009819795.1", taxa2 = "GCA_951394365.1",
  #            label = "GalloAnserformes", offset.text = 1, angle = 280+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_907165065.1", taxa2 = "GCF_036169615.1",
  #            label = "Columbea", offset.text = 1, angle = 296,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCF_017976375.1", taxa2 = "GCA_036417535.1",
  #            label = "Otidimorphae", offset.text = 1, angle = 295+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_009819775.1", taxa2 = "GCA_901699155.2",
  #            label = "Caprimulgimorphae", offset.text = 1, angle = 296+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_030867145.1", taxa2 = "GCA_964199755.1",
  #            label = "Core waterbirds", offset.text = 1, angle = 310+5,
  #            color = "#000") +
  # geom_strip(taxa1 = "GCA_031877795.1", taxa2 = "GCF_027579445.1",
  #            label = "Core landbirds", offset.text = 1, angle = 335+5,
  #            color = "#000") +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm"))

#UR_tree_bl
#ggsave("./UR_tree.pdf", plot = UR_tree, width = 15, height = 15)


# Uniqueness Ratio Tree
UR_tree = ggtree(tree, 
                 color = "black", 
                 size = 0, 
                 branch.length="none"
) %<+% UR_states +
  geom_tree(aes(color = state), continuous = 'colour', size=1) +
  scale_color_gradientn(colours=rev(c("red", 'orange', 'green', 'cyan', 'blue')),
                        breaks=c(0.2, 0.4, 0.6),
                        limits=c(0, 0.8))+
  layout_circular() + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(24)) +
  guides(color = guide_colorbar("Repeat\nRatio")) +
  geom_strip(taxa1 = "GCA_964187855.1", taxa2 = "GCA_964198595.1",
             label = "Agnatha", offset.text=1, angle = 0,
             color = "#000") +
  geom_strip(taxa1 = "GCA_035084275.1", taxa2 = "GCA_964213995.1",
             label = "Chondrichthyes", offset.text=1, angle = 7,
             color = "#000") +
  geom_strip(taxa1 = "GCF_900747795.2", taxa2 = "GCF_904848185.1",
             label = "Actinopterygii", offset.text=1, angle = 65,
             color = "#000") +
  geom_strip(taxa1 = "GCA_040939525.1", taxa2 = "GCF_037176945.1",
             label = "Sarcopterygii", offset.text=10, angle = 289.5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_901001135.1", taxa2 = "GCF_902459505.1",
             label = "Gymnophiona", offset.text=11, angle = 291,
             color = "#000") +
  geom_strip(taxa1 = "GCF_040938575.1", taxa2 = "GCA_964261635.1",
             label = "Caudata", offset.text=7, angle = 291.5+7,
             color = "#000") +
  geom_strip(taxa1 = "GCA_040206685.1", taxa2 = "GCF_040894005.1",
             label = "Anura", offset.text=5, angle = 300+5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_004115215.2", taxa2 = "GCF_015852505.1",
             label = "Monotremes", offset.text=9, angle = 310+5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_027887165.1", taxa2 = "GCF_011100635.1",
             label = "Marsupials", offset.text=8, angle = 315+5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_030445035.2", taxa2 = "GCA_023851605.1",
             label = "Xenartha", offset.text=7, angle = 319+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_043290085.1", taxa2 = "GCF_030014295.1",
             label = "Afrotheria", offset.text=7, angle = 320+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_026018925.1", taxa2 = "GCA_949316315.1",
             label = "Supraprimates", offset.text=11, angle = 340+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_964194135.1", taxa2 = "GCF_011762595.1",
             label = "Laurasiatheria", offset.text=11, angle = 30+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_046126795.1", taxa2 = "GCA_964234715.1",
             label = "Squamata", offset.text=8, angle = 65+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_045364815.1", taxa2 = "GCA_040894355.2",
             label = "Testudines (Pleurodira)", offset.text=16, angle = 79+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_033958435.1", taxa2 = "GCA_965140285.1",
             label = "Testudines (Cryptodira)", offset.text=16, angle = 87+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_030020295.1", taxa2 = "GCF_030867095.1",
             label = "Crocodylia", offset.text = 1, angle = 270+5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_040807025.1", taxa2 = "GCF_036370855.1",
             label = "Palaeognathae", offset.text = 1, angle = 275+5,
             color = "#000") +
  geom_strip(taxa1 = "GCF_009819795.1", taxa2 = "GCA_951394365.1",
             label = "GalloAnserformes", offset.text = 1, angle = 280+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_907165065.1", taxa2 = "GCF_036169615.1",
             label = "Columbea", offset.text = 1, angle = 296,
             color = "#000") +
  geom_strip(taxa1 = "GCF_017976375.1", taxa2 = "GCA_036417535.1",
             label = "Otidimorphae", offset.text = 1, angle = 295+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_009819775.1", taxa2 = "GCA_901699155.2",
             label = "Caprimulgimorphae", offset.text = 1, angle = 296+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_030867145.1", taxa2 = "GCA_964199755.1",
             label = "Core waterbirds", offset.text = 1, angle = 310+5,
             color = "#000") +
  geom_strip(taxa1 = "GCA_031877795.1", taxa2 = "GCF_027579445.1",
             label = "Core landbirds", offset.text = 1, angle = 335+5,
             color = "#000") +
  theme(plot.margin = unit(c(0,0,0,0), units = "cm"))

#UR_tree
#ggsave("./UR_tree.pdf", plot = UR_tree, width = 15, height = 15)

ggsave("./UR_nooutgroups.pdf", cowplot::plot_grid(UR_tree, UR_tree_bl), width = 30, height = 15)
