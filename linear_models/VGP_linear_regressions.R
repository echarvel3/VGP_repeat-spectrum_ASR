library(phytools)
library(dplyr)
library(data.table)

setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/")
source("./comparision_to_other_genomes/VGP_genome_comparison.R")
setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/linear_models/")

# =========================================== #
#          PHYLOGENETIC REGRESSION            #
# =========================================== #

# Formatting Data

existing_RepToLength = all_UR_states[,c("node", "ScientificName.x", "state.y", "Extended lineage")]
existing_RepToLength = merge(x = existing_RepToLength, all_GL_states[,c("node", "state.y")])[!is.na(ScientificName.x)]
existing_RepToLength = merge(x = existing_RepToLength, all_HCRM_states[,c("node", "state.y")])[!is.na(ScientificName.x)]
names(existing_RepToLength) = c("node", "scientific_name", "RR", "extended_lineage", "GL", "HCRM")
#existing_RepToLength$RR = 1 - existing_RepToLength$UR

vgp_RepToLength = all_UR_states[,c("node", "ScientificName.x", "state.x", "Extended lineage")]
vgp_RepToLength = merge(x = vgp_RepToLength, all_GL_states[,c("node", "state.x")])[!is.na(ScientificName.x)]
vgp_RepToLength = merge(x = vgp_RepToLength, all_HCRM_states[,c("node", "state.x")])[!is.na(ScientificName.x)]
names(vgp_RepToLength) = c("node", "scientific_name", "RR", "extended_lineage", "GL", "HCRM")
#vgp_RepToLength$RR = 1 - vgp_RepToLength$UR

# ------------------- PHYLOGENETIC INDEPENDENT CONTRASTS ----------------------- #
#Labeling internal nodes
tree$node.label = paste("N",1:length(tree$node.label),sep="")

# Running PIC
vgp_RR_pic = pic(x = vgp_UR_annots, phy = tree)
vgp_GL_pic = pic(x = vgp_GL_annots, phy = tree)
pic_table = data.table(RR_pic = vgp_RR_pic, GL_pic = vgp_GL_pic, dataset = "VGP Assemblies", node=names(vgp_GL_pic))

existing_RR_pic = pic(x = existing_UR_annots, phy = tree)
existing_GL_pic = pic(x = existing_GL_annots, phy = tree)
pic_table = rbind(pic_table, 
                  data.table(RR_pic = existing_RR_pic, GL_pic = existing_GL_pic, dataset = "Existing Assemblies", node=names(vgp_GL_pic))
                  )

#Counting number of leaves under each node
pic_table$leafcount=vapply(pic_table$node, function(x) extract.clade(tree, node = x)$Nnode+1 , FUN.VALUE = 1)

#Labeling nodes of interest
pic_table[node == "N165", nodelabel := "Bonobo - Chimpanzee"]
pic_table[node == "N5", nodelabel := "Amniotes"]
pic_table[node == "N4", nodelabel := "Osteichthyes"]

#Plot (with no ouliers)
pic_correlation_plot= 
  ggplot(pic_table %>%
           filter(GL_pic < 2.5e+10 ), # has an outlier
         aes(x = RR_pic, y = GL_pic,fill=leafcount)) +
  stat_smooth(method = "lm", se = 0,color="black") +
  geom_vline(xintercept = 0,linetype="dotted",color="grey40")+
  geom_hline(yintercept = 0,linetype="dotted",color="grey40")+
  geom_point(shape = 21,color="black") +
  ggrepel::geom_label_repel(aes(label = ifelse(node %in% c("N4", "N5", "N165"), 
                                               yes = nodelabel, no = NA),
                                color = ifelse(node %in% c("N4", "N5"),
                                               yes = "white", no = "black")), show.legend = F) +
  scale_color_manual(values = c("black", "white")) +
  facet_wrap(.~dataset, scales = "free_x") +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), 
                   method = "spearman", 
                   label.sep = "", 
                   cor.coef.name = 'rho', size = 5,color="black", position = position_dodge2()) + #only cherries p = 0.77 & p = 0.92
  theme_classic() + ylab("Genome Length PIC") + xlab("Repeated 31-mers PIC\n") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        title = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.1, color="grey80"),
        legend.position =c(.9,.19),
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=15),
        legend.key.spacing.x = unit(2.0, "cm"),
        strip.text = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_fill_viridis_b(breaks=c(2,3,5,20,50,200),direction = -1,name="# Leaves")

pic_correlation_plot

# Node labeling test figure...
# treeplot = ggtree(tree) %<+% all_HCRM_states +
#   geom_label(aes(label = label), 
#              size = 3) + 
#   geom_point2(aes(subset = label %in% c("N164", "N167", "N3", "N4")), 
#               color = "blue", 
#               size = 3, alpha = 0.5) + 
#   layout_circular() +  
#   geom_tippoint(aes(color = Superorder )) +
#   #geom_tiplab() +
#   geom_strip(taxa1 = "Myxine_glutinosa", taxa2 = "Petromyzon_marinus",
#              label = "Agnatha", #offset.text=1, #angle = 1, 
#              color = "#000") +
#   geom_strip(taxa1 = "Stegostoma_tigrinum", taxa2 = "Carcharodon_carcharias",
#              label = "Chondrichthyes", #offset.text=1, #angle = 1, 
#              color = "#000") +
#   geom_strip(taxa1 = "Acipenser_ruthenus", taxa2 = "Scatophagus_argus", 
#              label = "Actinopterygii", #offset.text=1, #angle = 60,
#              color = "#000") + 
#   geom_strip(taxa1 = "Protopterus_annectens", taxa2 = "Latimeria_chalumnae", 
#              label = "Sarcopterygii", #offset.text=8, #angle = 288,
#              color = "#000") +
#   geom_strip(taxa1 = "Ascaphus_truei", taxa2 = "Engystomops_pustulosus",
#              label = "Anura", #offset.text=5, #angle = 290,
#              color = "#000") +
#   geom_strip(taxa1 = "Ornithorhynchus_anatinus", taxa2 = "Tachyossus_aculeatus",
#              label = "Monotremes", #offset.text=9, #angle = 300,
#              color = "#000") +
#   geom_strip(taxa1 = "Monodelphis_domestica", taxa2 = "Sminthopsis_crassicaudata",
#              label = "Marsupials", #offset.text=8, #angle = 310,
#              color = "#000") +
#   geom_strip(taxa1 = "Dasypus_novemcinctus", taxa2 = "Choloepus_didactylus",
#              label = "Xenartha", #offset.text=7, #angle = 319,
#              color = "#000") +
#   geom_strip(taxa1 = "Dugong_dugon", taxa2 = "Elephas_maximus_indicus",
#              label = "Afrotheria", #offset.text=7, #angle = 320,
#              color = "#000") +
#   geom_strip(taxa1 = "Tupaia_tana", taxa2 = "Marmota_flaviventris",
#              label = "Supraprimates", #offset.text=11, #angle = 340,
#              color = "#000") +
#   geom_strip(taxa1 = "Erinaceus_europaeus", taxa2 = "Globicephala_melas",
#              label = "Laurasiatheria", #offset.text=11, #angle = 30,
#              color = "#000") +
#   geom_strip(taxa1 = "Tiliqua_scincoides", taxa2 = "Podarcis_muralis",
#              label = "Squamata", #offset.text=8, #angle = 65,
#              color = "#000") +
#   geom_strip(taxa1 = "Podocnemis_expansa", taxa2 = "Podocnemis_expansa",
#              label = "Testudines (Pleurodira)", #offset.text=16, #angle = 65,
#              color = "#000") +
#   geom_strip(taxa1 = "Carettochelys_insculpta", taxa2 = "Chelonia_mydas",
#              label = "Testudines (Cryptodira)", #offset.text=16, #angle = 70,
#              color = "#000") +
#   geom_strip(taxa1 = "Gavialis_gangeticus", taxa2 = "Alligator_mississippiensis",
#              label = "Crocodylia", #offset.text=9, #angle = 70,
#              color = "#000") +
#   geom_strip(taxa1 = "Eudromia_elegans", taxa2 = "Rhea_pennata",
#              label = "Palaeognathae", #offset.text = 10, #angle = 260 + 170,
#              color = "#000") +
#   geom_strip(taxa1 = "Gallus_gallus", taxa2 = "Anser_brachyrhynchus",
#              label = "GalloAnserformes", #offset.text = 10, #angle = 260+180,
#              color = "#000") +
#   geom_strip(taxa1 = "Cuculus_canorus", taxa2 = "Chlamydotis_macqueenii",
#              label = "Otidimorphae", #offset.text = 1, #angle = 270,
#              color = "#000") +
#   geom_strip(taxa1 = "Phoenicopterus_ruber", taxa2 = "Caloenas_nicobarica",
#              label = "Columbea", #offset.text = 1, #angle = 280,
#              color = "#000") +
#   geom_strip(taxa1 = "Caprimulgus_europaeus", taxa2 = "Podargus_strigoides",
#              label = "Caprimulgimorphae", #offset.text = 1, #angle = 285,
#              color = "#000") +
#   geom_strip(taxa1 = "Opisthocomus_hoazin", taxa2 = "Larus_argentatus",
#              label = "Core waterbirds", #offset.text = 1, #angle = 305,
#              color = "#000") +
#   geom_strip(taxa1 = "Strix_aluco", taxa2 = "Zonotrichia_albicollis",
#              label = "Core landbirds", #offset.text = 1, #angle = 335,
#              color = "#000") +
#   theme(plot.margin = unit(c(0,0,0,0), units = "cm"))
# 
# ggsave(filename = "test.pdf", plot = treeplot, width = 40, height = 40)


# ------------------- PHYLOGENETIC GENERALIZED LEAST SQUARES ----------------------- #

tree = read.tree("../roadies_v1.1.16b.nwk") 
#Read Taxonomy Annotations
merged_annots = fread("../annotations_vgp.txt")
#Read VGP Respect Results
respect_results = fread("../full_VGP_analyses/respect_full-VGP_parameters.txt")
respect_results$sample = sub(".hist", "", respect_results$sample)
#Format Data
merged_annots = merge(x = respect_results,  y = merged_annots, by.x = "sample", by.y = "Assc.",)
vgp_RepToLength = merge(merged_annots, extend_lin, by.x = "Family", by.y = "Family Scientific Name")
vgp_RepToLength = vgp_RepToLength[,c("sample", "genome_length", "uniqueness_ratio", "Extended lineage","Family","ScientificName", "Superorder")]
vgp_RepToLength$RR = 1 - vgp_RepToLength$uniqueness_ratio
names(vgp_RepToLength) = c("sample", "genome_length", "UR", "extended_lineage", "family","scientific_name", "superorder", "RR")
vgp_RepToLength$GL = vgp_RepToLength$genome_length / 1000000
vgp_RepToLength$RRb = vgp_RepToLength$RR
vgp_RepToLength$RR = vgp_RepToLength$RRb +runif(min = 0,max= 0.01,n=nrow(vgp_RepToLength))

# Make linear model given lineage label
generate_lm = function(x){
  tips =  vgp_RepToLength[`extended_lineage` != x, "sample"]
  small_tree = drop.tip(phy = tree, tip = as.character(c(tips$sample)))
  species = small_tree$tip.label
  V = ape::corBrownian(phy = small_tree, form = ~sample)
  nlme::gls(GL ~ RR, correlation = V, data = vgp_RepToLength[sample %in% species])
}

#Make Linear model for each Lineage
Birds_lm = generate_lm("Birds")
Turtles_lm = generate_lm("Turtles")
Lepidosauria_lm = generate_lm("Lepidosauria")
Mammals_lm = generate_lm("Mammals")
Amphibians_lm = generate_lm("Amphibians")
Ray_lm = generate_lm("Ray-finned fishes")
Cartilaginous_lm = generate_lm("Cartilaginous fishes")
Cyclostomes_lm = generate_lm("Cyclostomes")

birds = data.table(Birds_lm$fitted, 
                   "Birds", 
                   vgp_RepToLength[extended_lineage == "Birds", c("superorder", "family", "scientific_name", "RR", "GL")])
turtles = data.table(Turtles_lm$fitted, 
                     "Turtles", 
                     vgp_RepToLength[extended_lineage == "Turtles",  c("superorder", "family", "scientific_name", "RR", "GL")])
lepidosauria = data.table(Lepidosauria_lm$fitted, 
                          "Lepidosauria", 
                          vgp_RepToLength[extended_lineage == "Lepidosauria",  c("superorder", "family", "scientific_name", "RR", "GL")])
mammals = data.table(Mammals_lm$fitted, 
                     "Mammals", 
                     vgp_RepToLength[extended_lineage == "Mammals",  c("superorder", "family", "scientific_name", "RR", "GL")])
amphibians = data.table(Amphibians_lm$fitted, 
                        "Amphibians", 
                        vgp_RepToLength[extended_lineage == "Amphibians",  c("superorder", "family", "scientific_name", "RR", "GL")])
ray_finned = data.table(Ray_lm$fitted, 
                        "Ray-finned fishes", 
                        vgp_RepToLength[extended_lineage == "Ray-finned fishes",  
                                        c("superorder", "family", "scientific_name", "RR", "GL")])
cartilaginous = data.table(Cartilaginous_lm$fitted,
                           "Cartilaginous fishes", 
                           vgp_RepToLength[extended_lineage == "Cartilaginous fishes",  c("superorder", "family", "scientific_name", "RR", "GL")])
cyclostomes = data.table(Cyclostomes_lm$fitted, 
                         "Cyclostomes", 
                         vgp_RepToLength[extended_lineage == "Cyclostomes",  c("superorder", "family", "scientific_name", "RR", "GL")])


#Re-run regression with no outliers

outliers = c(mammals[(GL)/V1 > 1.5, scientific_name], amphibians[(GL)/V1 > 1.5, scientific_name], ray_finned[(GL)/V1 > 1.5, scientific_name])
temp = vgp_RepToLength
vgp_RepToLength = vgp_RepToLength[!scientific_name %in% outliers]

Mammals_lm_noOutlier = generate_lm("Mammals")
Amphibians_lm_noOutlier = generate_lm("Amphibians")
Ray_lm_noOutlier = generate_lm("Ray-finned fishes")
vgp_RepToLength = temp

mammals = data.table(as.numeric(predict(Mammals_lm_noOutlier, data.frame(RR = vgp_RepToLength[extended_lineage == "Mammals", RR]))), 
                     "Mammals",
                     vgp_RepToLength[extended_lineage == "Mammals",  c("superorder", "family", "scientific_name", "RR", "GL")])

amphibians = data.table(as.numeric(predict(Amphibians_lm_noOutlier, data.frame(RR = vgp_RepToLength[extended_lineage == "Amphibians", RR]))),
                        "Amphibians", 
                        vgp_RepToLength[extended_lineage == "Amphibians",  c("superorder", "family", "scientific_name", "RR", "GL")])
ray_finned = data.table(as.numeric(predict(Ray_lm_noOutlier, data.frame(RR = vgp_RepToLength[extended_lineage == "Ray-finned fishes", RR]))), 
                        "Ray-finned fishes", 
                        vgp_RepToLength[extended_lineage == "Ray-finned fishes",  c("superorder", "family", "scientific_name", "RR", "GL")])


data = do.call("rbind", list(birds, turtles, lepidosauria, mammals, ray_finned, amphibians, cyclostomes, cartilaginous))
data$V4 = factor(data$V2,labels=c("Amphibians","Amniotes","Fishes","Fishes","Amniotes","Amniotes","Fishes","Amniotes"))

#################################

# LM Arrow Plot

lm_predictions = ggplot(data, aes(x = RR, y = GL, yend=V1, color = V2)) + 
  geom_segment(arrow = arrow(length = unit(2.5,"pt"),type = "open"),
               , alpha= 0.9)+
  facet_wrap(.~V4  , scales = "free", 
             nrow = 1) +
  ggrepel::geom_label_repel(aes(label = ifelse(scientific_name %in% outliers, scientific_name, NA)),
                            size=4, show.legend = F) +
  ylab("Genome Length") + xlab("Repeated 31-mer Ratio") +
  theme_classic() +
  guides(color = guide_legend("", override.aes = list(linewidth = 2, 
                                                      arrow = NULL))) +
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
  ),name="") +  scale_x_continuous(expand = expansion(mult = c(0.15, 0.1))) 

lm_predictions

ggsave(filename = "../figures/linear_modeling_results.pdf", 
       plot = cowplot::plot_grid(pic_correlation_plot, 
                                 lm_predictions,  
                                 ncol = 1, 
                                 labels = c("A", "B"), label_size = 30, label_y = 1.025,
                                 rel_heights = c(1,1.5)),
       height = 13, width = 13)
