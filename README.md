# VGP Repeat Spectrum Ancestral State Reconstruction  

This repository contains analyses of repeat spectrum ancestral state reconstruction (ASR) across genomes from the [Vertebrate Genomes Project (VGP)](https://vertebrategenomesproject.org/).  

---

## Contents  

### 1. `full_VGP_analyses/`  
ASR of all **568 genomes** in VGP Phase 1 using phylotools. Parameters inferred from using [RESPECT](https://github.com/shahab-sarmashghi/RESPECT) (k-mer size = 31).  

- `full_VGP_analyses/respect_full-VGP_parameters.txt`:

  Table of genome accessions with corresponding:  
  - **Uniqueness Ratio** (UR)  -- _number of unique k-mers / genome length_
  - **High Copy Repeats per Million Base Pairs** (HCRM) -- _The count of the 10 most repetitive k-mers, divided by genome length_
  - **Genome Length** (GL)

- `full_VGP_analyses/respect_full-VGP_spectra.txt`
  
  Table of genome accessions with *k*-mer counts for the first 50 *k*-mer frequencies.  

- `full_VGP_analyses/ASR/*.nwk`
  
  Directory contains phylogenetic tree with internal nodes containing values from a maximum likelihood ancestral state reconstructions (using a Brownian motion model with the _phylotools_ R package) for three parameters:  
  - RR (Repeat Ratio) [`RR = 1 - UR`] -- _full_VGP_UR-ASR.nwk_
  - HCRM -- _full_VGP_HCRM-ASR.nwk_
  - GL -- _full_VGP_GL-ASR.nwk_
    
- `full_VGP_analyses/VGP_ancestral_recon.R`
  - R script used for full VGP dataset ASR for all three parameters and corresponding figures.

- **Corresponding Figures:**  
  - ASR Trees for all RR, HCRM, and GL: `figures/VGP_full_ASR.pdf`  

---

### 2. `comparision_to_other_genomes/`  
Analysis of a **subset of 248 VGP genomes** with corresponding existing genome assemblies, enabling direct comparison between VGP assemblies and prior datasets.  
- `comparison_to_other_genomes/roadies_v1.1.4.subsampled.species-labels.nwk`
   - Phylogenetic tree constructed with full VGP data, but pruned to only keep assemblies with both VGP and existing assemblies.

- `comparision_to_other_genomes/respect_existing-VGP_parameters.txt`:

  Table of genome accessions with corresponding: UR, HCRM, and GL for **all existing assemblies** of VGP species.

- `comparision_to_other_genomes/respect_full-existing_spectra.txt`
  
  Table of genome accessions with *k*-mer counts for the first 50 *k*-mer frequencies for **all existing assemblies** of VGP species.
    
- `comparision_to_other_genomes/parameter_comparison.tsv`

  contains the following data for subset of assemblies being used for comparison:
    - Lineage
    - Scientific Name
    - Superorder
    - Genome accession id *(for both existing and VGP genomes)*
    - UR *(for both existing and VGP genomes)*
    - HCRM *(for both existing and VGP genomes)*
    - GL *(for both existing and VGP genomes)*
    - RR: *(for both existing and VGP genomes)*
    - Ratio of `(existing_RR)/(VGP_RR)`

- `comparison_to_other_genomes/ASR/*nwk`
  
  Directory contains phylogenetic tree with internal nodes containing values from a maximum likelihood ancestral state reconstructions.
  - `subset_VGP*nwk`: set of ASR trees with subsampled tree corresponding to VGP genome parameters
  - `subset_existing*nwk`: set of ASR trees with subsampled tree corresponding to existing genome parameters.`
 
- `comparison_to_other_genomes/VGP_genome_comparison.R`
  - R script used for running ASR for VGP and existing subset genomes, comparison scatter plots, and comparison trees.

- `comparison_to_other_genomes/VGP_spectra_comparison.R`
 - R Script for creating arrow plot comparing individual k-mer frequencies.

- **Corresponding Figures:**  
  - Scatter plots comparing parameter values: `figures/all_plots.pdf`
  - ASR for corresponding datasets using subsampled tree: `figures/*_subset_ASR.pdf`
  - Trees with % change in parameters across corresponding ASR: `figures/ASR_comparison_trees.pdf`
  - Tree with % change in RR with blue-red color scheme `figures/RR_comparison_tree.pdf`
  - Arrow plots showing change in k-mer counts (`figures/arrow_plot_kmer_count.pdf`) and portion of the genome represented by k-mer frequency (`figures/arrow_plot.pdf`).

### 3. `linear_models/`  
- `linear_models/VGP_linear_regressions.R`
 - R Script for running Phylogenetic Independent Contrasts (PIC) and phylogenetic regressions.

- **Corresponding Figures:**  
  - Panel A shows Spearman's Correlation of GL and RR PICs, panel B shows GL ~ RR phylogenetic regression results, and illustrates change in predicted and true GL values: `figures/linear_modeling_results.pdf`
