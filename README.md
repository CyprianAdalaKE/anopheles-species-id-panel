# Anopheles species-ID and Plasmodium detection amplicon panel

Analysis code for:

> **A targeted nanopore amplicon sequencing panel for simultaneous identification of six *Anopheles* species and *Plasmodium* detection.**
> Guya CA *et al.* (manuscript in preparation)

The panel is an 11-target nanopore multiplex: 10 nuclear loci for species identification across
*An. gambiae*, *An. coluzzii*, *An. arabiensis*, *An. funestus*, *An. stephensi* and *An. coustani*,
plus one *Plasmodium* 18S rDNA marker (Plasmo18S).

Raw sequence data are available from the NCBI Sequence Read Archive under BioProject
**PRJNA1517867**. This repository holds scripts, small input and metadata files, and summary
result tables only. It does not include reads, alignments or reference genomes.

## Setup

```bash
git clone https://github.com/CyprianAdalaKE/anopheles-species-id-panel.git
cd anopheles-species-id-panel
conda env create -f envs/anopheles-panel.yml   # minimap2, samtools, bcftools, Python stack
conda env create -f envs/kraken2-bracken.yml   # Kraken2, Bracken, seqtk
conda env create -f envs/phylo_env.yml         # IQ-TREE
conda env create -f envs/medaka.yml            # Medaka (host-18S polishing only)
```

The environment files were exported from the analysis machine at the time of release.

Scripts expect the per-species FASTQ folders (`gam_final/`, `colu_final/`, `arabiensis/`,
`fun_final/`, `coustani_combined/`, `steph_combined/`, `infectivity/`) and the AgamP4 reference
(`VectorBase-67_AgambiaePEST_Genome.fasta`) in the repository root. Point the scripts elsewhere
with environment variables:

| Variable | Used by | Meaning |
|---|---|---|
| `PANEL_DIR` | root scripts, `build_manifest.py` | Data root (default: repository root) |
| `KRAKEN_DB` | `run_kraken_bracken.sh`, `build_database.sh`, `run_downsample_experiment.sh` | Kraken2/Bracken database directory |
| `BRACKEN_BIN` | `run_kraken_bracken.sh`, `run_downsample_experiment.sh` | Bracken executable (default: `bracken`) |
| `PF_REF` | `align_infectivity.sh` | PlasmoDB-67 *P. falciparum* 3D7 genome FASTA |
| `PLASMODB_GENOMES` | `plasmo_validation/build_indexes.sh` | Directory with the four PlasmoDB-67 genome FASTAs |

## Workflow

### 1. Custom Kraken2/Bracken database: `kraken2_database/`

| Script | Purpose |
|---|---|
| `add_kraken_to_fasta.py` | Tags VectorBase/PlasmoDB FASTA headers with `kraken:taxid\|<taxid>` (Supplementary File 1) |
| `build_database.sh` | Adds the 6 *Anopheles* + 4 *Plasmodium* genomes (release 67) to the library, builds Kraken2 (k = 35, ℓ = 31) and Bracken (100 bp reads) |

### 2. Alignment and coverage (Figures 7–8, Supplementary Table S1, Supplementary Figure S2)

| Script | Purpose |
|---|---|
| `align_all_samples.sh` | minimap2 alignment of every sample to AgamP4 |
| `align_infectivity.sh` | Alignment of infectivity-experiment samples to *P. falciparum* 3D7 |
| `01_compute_coverage_v2.sh` | Per-locus depth over `nomadsSID.amplicons.bed` |
| `02_aggregate_and_plot_v2.py` | Aggregates coverage and draws the depth heatmap |
| `figures_v2/replot_depth_and_sens.py` | Log-scale depth boxplots (Figures 7–8) and sensitivity/specificity with Wilson 95% CIs (Figure 10, Table 4). Output filenames use the draft figure numbers (8, 9, 12) |

### 3. Taxonomic classification

| Script | Purpose |
|---|---|
| `run_kraken_bracken.sh` | Kraken2 + Bracken species-level classification of all samples |

### 3b. Species-identification accuracy: `species_id_accuracy/` (Figure 9, Table 4, Figure 10)

`confusion_matrix_sens_spec.py` builds the confusion matrix of panel calls against reference
species identification and computes one-vs-rest sensitivity and specificity per species with
Wilson 95% CIs, plus overall accuracy (445/478 = 93.1%). The input CSV has one row per sample
with columns `panel` and `reference_species_id`:

```bash
python species_id_accuracy/confusion_matrix_sens_spec.py species_id_accuracy/species_id_joined_k2_Bracken_updated_Jul26_latest.csv species_id_accuracy/out
```

`crosstab_original_colab.ipynb` is the original Google Colab notebook used for the manuscript.

### 4. Genotype-based validation: `pca_validation/` (Figures 11, 12, 14)

Run in this order:

1. `build_manifest.py`: builds the sample → species → BAM table (`manifest.csv`).
2. `call_variants.sh`: joint diploid SNP calling with bcftools at the 10 target loci.
3. `filter_and_matrix.sh`: filters to biallelic SNPs and writes the genotype matrix.
4. `pca.py` / `pca2.py`: PCA, for all species and for the *An. gambiae* complex only.
5. `locus_discriminatory_power.py`: per-locus Hudson Fst and cross-validated LDA.

`check_missingness.py`, `check_loadings.py` and `check_tail.py` are QC helpers. The missingness
check identifies the six high-missingness *An. coluzzii* samples excluded from the PCA.

### 5. Phylogeny: `reviewer_prep/phylo_tree/` (Figure 13)

`build_alignment.py` builds the SNP alignment, `run_iqtree.sh` runs IQ-TREE (GTR+G, 1,000
ultrafast bootstraps) and `render_tree.py` draws the tree. The final tree files
(`snp_tree.treefile`, `snp_tree.contree`, `snp_tree.iqtree`) are included.

### 6. Diagnostic SNPs: `reviewer_prep/diagnostic_snps/` (Table 5)

`find_diagnostic_snps.py` reports the strongest allele-frequency-differentiating site per locus.
Output: `diagnostic_snps.tsv`.

### 7. Plasmo18S alignment-based reclassification: `plasmo_validation/` (Table 7)

| Script | Purpose |
|---|---|
| `build_indexes.sh` | minimap2 indexes for the host, Plasmodium and host-18S references |
| `align_all.sh` | Naive competition: AgamP4 vs the combined 4-species Plasmodium reference |
| `realign_host.sh`, `align_18s_standalone.sh` | Host reference supplemented with public *An. gambiae* 18S sequences |
| `align_18s_empirical.sh`, `align_18s_polished.sh` | Empirical host-18S reference, before and after Medaka polishing |
| `classify_reads.py` | Per-read competitive classification (identity ≥ 0.85, coverage ≥ 0.70, margin 0.03) |

Host-18S reference files:

- `angambiae_18s_*.fasta`: public *An. gambiae* 18S sequences.
- `draft_host18s.fasta`, `host18s_empirical_r1.fasta`, `host18s_empirical.fasta`: iterative
  consensus from mosquito-only negative-control reads.
- `host18s_polished.fasta`: the empirical reference polished with Medaka 2.2.0.

These references were built interactively and not as a scripted step, so the final FASTAs are
provided directly.

### 8. Depth vs accuracy downsampling: `reviewer_prep/depth_accuracy/` (Figure 19, Table 9)

`select_samples.py` picks 48 samples, 8 per species and spanning each species' native depth
range. `run_downsample_experiment.sh` downsamples each one with seqtk (seed 100) to 50–5,000
reads and reclassifies it with Kraken2/Bracken. Results are in `depth_experiment_results.tsv`.

## Resources

- `nomadsSID.amplicons.bed`: target amplicon coordinates on AgamP4.
- `resources/primers.tsv`: primer sequences for the 10 species-ID loci (Table 2).

## License

MIT. See [LICENSE](LICENSE).
