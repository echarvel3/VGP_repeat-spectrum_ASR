# VGP Repeat Spectrum Ancestral State Reconstruction  

This repository contains analyses of repeat spectrum ancestral state reconstruction (ASR) across genomes from the [Vertebrate Genomes Project (VGP)](https://vertebrategenomesproject.org/).  

---

## Contents  

### 1. `full_VGP_analyses/`  
Ancestral state reconstruction (ASR) of all **568 genomes** in the VGP using [RESPECT](https://github.com/shahab-sarmashghi/RESPECT) (k-mer size = 31).  

- **respect_full-VGP_parameters.txt**  
  Table of genome accessions with:  
  - Uniqueness ratios (UR)  
  - HCRM (High Copy Repeats per Million)  
  - Genome lengths  

- **respect_full-VGP_spectra.txt**  
  Table of genome accessions with *k*-mer counts for the first 50 *k*-mer frequency bins.  

- **ASR/**  
  Maximum likelihood ancestral state reconstructions (using a Brownian motion model with the _phylotools_ R package) for three parameters:  
  - UR (Uniqueness Ratio; plotted as `RR = 1 - UR`)  
  - HCRM (High Copy Repeats per Million)
  - GL (Genome Length)  

- **Figures**  
  - Annotated phylogenetic tree: `figures/full_VGP_ASR.pdf`  
  - Panel **B** of `figures/kmerrepeat.pdf`  

---

### 2. `comparision_to_other_genomes/`  
Analysis of a **subset of 257 VGP genomes** with corresponding existing whole-genome assemblies, enabling direct comparison between VGP assemblies and prior datasets.  

- **Data**  
  - For each species in the analysis `comparision_to_other_genomes/parameter_comparison.tsv` contains:
    - Lineage
    - Scientific Name
    - Superorder
    - Genome accession id *(for both existing and VGP genomes)*
    - UR *(for both existing and VGP genomes)*
    - HCRM *(for both existing and VGP genomes)*
    - GL *(for both existing and VGP genomes)*
    - RR: *(for both existing and VGP genomes)*
    - Ratio of `(existing_RR)/(VGP_RR)`

- **ASR**  
  - Performed for both existing and VGP genomes using the same phylogeny, with tip values set for each set of assemblies.
  - Annotated phylogenies: `comparision_to_other_genomes/ASR/`  

- **Figures**  
  - Comparison of parameter distributions: `/figures/all_plots.pdf`
  - ASR for existing assemblies: `figures/subset_existing_ASR.pdf`  
  - Panel **A** and **C** of `figures/kmerrepeat.pdf`  
