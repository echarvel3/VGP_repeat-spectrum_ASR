
library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggtree)
library(dplyr)

# Helper function:
split_row <- function(x){
  z = paste(collapse = "_", strsplit(x, split = "_")[[1]][1:2])
}

# =========================================== #
#              Data Pre-processing            #
# =========================================== #

setwd("~/Desktop/VGP_analyses/VGP_repeat-spectrum_ASR/comparision_to_other_genomes/")
# Read Subsampled Tree (have both VGP and existing assemblies)
tree = read.tree("../roadies_v1.1.16b.nwk") 
# tree = drop.tip(tree, tip = c("Sminthopsis_crassicaudata","Shinisaurus_crocodilurus","Scardinius_erythrophthalmus","Tupaia_tana","Macrotis_lagotis","Aulostomus_maculatus","Salvelinus_alpinus","Spinachia_spinachia","Oenanthe_melanoleuca"))
# tree$node.label = paste("N",1:length(tree$node.label),sep="")

# ------------ Reading VGP Data ------------- #
# Reading VGP annotations
vgp_annots = fread("../annotations_vgp.txt")
extend_lin = unique.data.frame(fread("../VGP_extended_lin.tsv")[,c("Family Scientific Name", "Extended lineage")])
vgp_annots = merge(vgp_annots, extend_lin, by.x = "Family", by.y = "Family Scientific Name", all.x = T, all.y = F, )
vgp_annots$ScientificName = gsub(pattern = " ", replacement = "_", x =vgp_annots$ScientificName)

# Reading VGP Assemblies Respect results
vgp_parameters = fread("../full_VGP_analyses/respect_full-VGP_parameters.txt")
vgp_parameters$accession = gsub(pattern = ".hist", replacement = "", x = vgp_parameters$sample)

vgp_merged = as.data.table(merge(vgp_parameters, vgp_annots, by.x = 'accession', by.y = 'Assc.'))

# Merge all VGP data
#vgp_merged = as.data.table(subset(x = vgp_merged, subset = vgp_merged$ScientificName %in% tree$tip.label))
vgp_merged$genome_length = as.numeric(vgp_merged$genome_length)
vgp_merged$repeat_ratio = 1 - as.numeric(vgp_merged$uniqueness_ratio)
#vgp_merged$GL_mbp = as.numeric(vgp_merged$genome_length)/1000000

# Change to phylogeny tips from genome accession to species name
label_map = data.frame(old_label = vgp_annots$`Assc.`, new_label = vgp_annots$ScientificName)
tree$tip.label[match(label_map$old_label, tree$tip.label)] <- label_map$new_label

# ------- Reading 'Existing Genome' Data -------- #

# Reading Existing Assembly annotations
existing_annots = read.csv("./curated_low_qual_species.tsv")#[,c("Assc.", "ScientificName")]

# Reading Existing Assembly Respect results
existing_parameters = fread("./respect_full-existing_parameters.txt")
existing_parameters$accession= sapply(X = existing_parameters$sample, FUN = split_row)
existing_annots$ScientificName = existing_annots$scientific_name

# Merge all Existing Assembly data
existing_merged = as.data.table(merge(existing_parameters, existing_annots, by.x = 'accession', by.y = 'assembly_id'))
existing_merged[ScientificName == "Gorilla gorilla gorilla", ScientificName := "Gorilla gorilla"]
existing_merged[ScientificName == "Phoenicopterus ruber ruber", ScientificName := "Phoenicopterus ruber"]
existing_merged[ScientificName == "Canis lupus", ScientificName := "remove"]
existing_merged[ScientificName == "Oenanthe melanoleuca", ScientificName := "remove"]
existing_merged[ScientificName == "Rhinophrynus dorsalis", ScientificName := "remove"]
existing_merged[ScientificName == "Alca torda", ScientificName := "Alca Torda"]
existing_merged[ScientificName == "Mustela nivalis", ScientificName := "Mustela nivalis vulgaris"]
existing_merged[ScientificName == "Monodon monoceros", ScientificName := "Monodon monocero"]
existing_merged[ScientificName == "Lonchura striata", ScientificName := "Lonchura striata domestica"]
existing_merged[ScientificName == "Guaruba guarouba", ScientificName := "Guaruba guaruba"]
existing_merged[ScientificName == "Salvelinus alpinus oquassa", ScientificName := "Salvelinus alpinus"]
existing_merged[ScientificName == "Rhamphochromis sp. 'chilingali'", ScientificName := "Rhamphochromis sp. chilingali"]
existing_merged[ScientificName == "Lycocorax pyrrhopterus obiensis", ScientificName := "Lycocorax pyrrhopterus"]
existing_merged = existing_merged[!accession %in% c("GCA_019141155.1", 
                                                    "GCA_021292165.1",
                                                    "GCA_026018925.1",
                                                    "GCA_036784965.1",
                                                    "GCA_037893015.1",
                                                    "GCA_048126635.1",
                                                    "GCA_048301465.1",
                                                    "GCA_048593235.1",
                                                    "GCF_029582105.1",
                                                    "GCA_024453875.1"
)]
existing_merged$genome_length = as.numeric(existing_merged$genome_length)
existing_merged$repeat_ratio = 1 - as.numeric(existing_merged$uniqueness_ratio)
#existing_merged$GL_mbp = as.numeric(existing_merged$genome_length)/1000000
existing_merged$ScientificName = gsub(pattern = " ", replacement = "_", x = existing_merged$ScientificName)
existing_merged = as.data.table(subset(x = existing_merged, subset = existing_merged$ScientificName %in% tree$tip.label))

tree = drop.tip(tree, tree$tip.label[!tree$tip.label %in% existing_merged$ScientificName])
vgp_merged = vgp_merged[ScientificName %in% tree$tip.label]
existing_merged = existing_merged[ScientificName %in% tree$tip.label]
