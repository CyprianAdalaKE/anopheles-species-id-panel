# Anopheles species-ID and Plasmodium detection amplicon panel

Analysis code for:

> **A targeted nanopore amplicon sequencing panel for identification of six *Anopheles* species, with an evaluation of integrated *Plasmodium* detection.**
> Guya CA, Odhiambo CA, Odhiambo BM, Omoke D, Millanoi S, Polo B, Mwenda M, Awandu S, Onguru D, Bridges DJ, Ochomo E, Hendry J. Submitted to *BMC Genomics* (2026).

The panel is an 11-target nanopore multiplex: 10 nuclear loci for species identification across
*An. gambiae*, *An. coluzzii*, *An. arabiensis*, *An. funestus*, *An. stephensi* and *An. coustani*,
plus one *Plasmodium* 18S rDNA marker (Plasmo18S).

Raw sequence data are available from the NCBI Sequence Read Archive under BioProject
**PRJNA1517867**. This repository holds scripts, small input and metadata files, and summary
result tables only. It does not include reads, alignments or reference genomes.

## Manuscript figures and tables

All results in the submitted manuscript use the evaluation set of 349 specimens
(72 *An. arabiensis*, 76 *An. gambiae*, 77 *An. funestus*, 75 *An. stephensi*, 21 *An. coluzzii*,
28 *An. coustani*). Every read was checked against the original basecalled output of its run. Of
374 specimens sequenced, 25 were excluded by the quality-control rule below (see `qc_checks/`).

| Item | Script | Output / data |
|---|---|---|
| Fig. 1 (study-site map) | `figure_scripts/fig1_map/fig1_map.py` | `figure1_study_sites.png` |
| Fig. 2 (workflow) | `figure_scripts/fig2_overview.py` | |
| Fig. 3 (ortholog conservation screen) | `figure_scripts/fig3.py` | |
| Table 3, Fig. 4 (read yield, read length) | `figure_scripts/fig6_table3_yield_v2.py` | |
| Figs 5–6 (per-locus depth) | `figure_scripts/fig7_8_depth.py` | |
| Fig. 7, Table 4 (confusion matrix, sensitivity/specificity; 88.0%, 307/349) | `verified_set/run_k2_verified6.sh`, `verified_set/accuracy6.py`, `figure_scripts/fig7_confusion.py` | `verified_set/per_sample_calls.tsv`, `verified_set/accuracy_summary.txt` |
| Figs 8–10 (PCA, Fst), LDA | `pca_validation/verified_set/` (`run_genotyping.sh`, `pca.py`, `pca2.py`, `locus_discriminatory_power.py`, `lda_complex.py`), `figure_scripts/fig11_12_14.py` | `pca_validation/verified_set/results/` |
| Fig. 11 (ML tree; IQ-TREE, up to 5,000 iterations) | `reviewer_prep/phylo_tree/` (`build_alignment.py`, `run_iqtree.sh`, `fig13_tree_nm5000.py`, `tree_pure_clades.py`) | `snp_tree_nm5000.*`, `clade_summary_nm5000.tsv` |
| Table 5 (diagnostic SNPs) | `pca_validation/verified_set/find_diagnostic_snps.py` | |
| Tables 6–7, Fig. 12 (Plasmo18S) | `plasmo_validation/`, `figure_scripts/plasmo_combined.py` | |
| Table 8, Fig. 13 (downsampling) | `reviewer_prep/depth_accuracy/`, `figure_scripts/fig19_table9_depth.py` | |
| Table S1 (primer-site mismatches, interspecies differences) | `figure_scripts/tableS1_primer_site_divergence.py` | |
| Table S2, Fig. S4 (samples below 500×) | `figure_scripts/figS2_tableS1.py` | |
| Tables S4–S5 (read-level identity, per-locus accuracy, Bracken read length) | `reviewer_analyses/` (`make_tables_S4_S5.py`) | `reviewer_analyses/summary.txt` and `*.tsv` |
| Quality control, carry-over, per-run accuracy, read-identity (novelty) rule; Additional file 2 | `qc_checks/` | `qc_checks/*.tsv`, `manuscript_numbers_349.txt` |
| Table 3 / Table 4 CIs, Table S2 | `figure_scripts/` | `figure_scripts/tables/` |
| SRA submission metadata | `sra_submission/` | `files_manifest.tsv`, `md5_source_files.txt` |

Scripts carry the absolute paths of the analysis machine in variables at the top of each file;
edit these, or the environment variables below, to rerun them elsewhere. Folders and scripts not
listed above (for example `species_id_accuracy/`, `figures_v2/`, `pca_validation/*.py` at the top
level) are from earlier analyses of a larger, unverified sample set and are kept for the record;
they do not produce the numbers in the submitted manuscript.

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

### 2. Alignment and coverage (Figs 5–6, Table S2, Fig. S4)

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

### 3b. Species-identification accuracy

The manuscript numbers (Table 4, Fig. 7) come from `verified_set/` (see the table above). `species_id_accuracy/` below is an earlier analysis of a larger, unverified sample set and is superseded.

`confusion_matrix_sens_spec.py` builds the confusion matrix of panel calls against reference
species identification and computes one-vs-rest sensitivity and specificity per species with
Wilson 95% CIs, plus overall accuracy (445/478 = 93.1%). The input CSV has one row per sample
with columns `panel` and `reference_species_id`:

```bash
python species_id_accuracy/confusion_matrix_sens_spec.py species_id_accuracy/species_id_joined_k2_Bracken_updated_Jul26_latest.csv species_id_accuracy/out
```

`crosstab_original_colab.ipynb` is the original Google Colab notebook used for the manuscript.

### 4. Genotype-based validation: `pca_validation/` (Figs 8–10)

The manuscript uses `pca_validation/verified_set/` (349 samples, 988 SNPs; `run_genotyping.sh` replaces steps 1–3). The top-level scripts are the original versions.

Run in this order:

1. `build_manifest.py`: builds the sample → species → BAM table (`manifest.csv`).
2. `call_variants.sh`: joint diploid SNP calling with bcftools at the 10 target loci.
3. `filter_and_matrix.sh`: filters to biallelic SNPs and writes the genotype matrix.
4. `pca.py` / `pca2.py`: PCA, for all species and for the *An. gambiae* complex only.
5. `locus_discriminatory_power.py`: per-locus Hudson Fst and cross-validated LDA.

`check_missingness.py`, `check_loadings.py` and `check_tail.py` are QC helpers. In the 349-sample
set no sample exceeds 50% missing genotype calls, so all are retained.

### 5. Phylogeny: `reviewer_prep/phylo_tree/` (Fig. 11)

`build_alignment.py` builds the SNP alignment and `run_iqtree.sh` runs IQ-TREE (GTR+G, 1,000
ultrafast bootstraps). The manuscript tree was run with up to 5,000 iterations
(`iqtree -s snp_alignment.fasta -st DNA -m GTR+G -B 1000 -nm 5000 -seed 669192 -T 4 -pre snp_tree_nm5000`);
on the 349-sample, 988-SNP alignment. Ultrafast bootstrap did not reach the default convergence criterion
(split-frequency correlation 0.86), so support values are approximate. `fig13_tree_nm5000.py` draws the tree
and writes `clade_summary_nm5000.tsv`; `tree_pure_clades.py` reports the largest species-pure clade per species
(*An. stephensi* 71/75 and *An. funestus* 50/77, both 100% support).

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

### 8. Depth vs accuracy downsampling: `reviewer_prep/depth_accuracy/` (Fig. 13, Table 8)

`select_samples.py` picks 48 samples, 8 per species and spanning each species' native depth
range. `run_downsample_experiment.sh` downsamples each one with seqtk (seed 100) to 50–5,000
reads and reclassifies it with Kraken2/Bracken. Results are in `depth_experiment_results.tsv`.

### 9. Read-level and locus-level analyses: `reviewer_analyses/` (Tables S4–S5)

| Script | Purpose |
|---|---|
| `extract_amplicons.py` | Amplicon sequence of each species from the ortholog MSAs |
| `map_msa_amplicons.sh` | Locates those amplicons in each species' genome (`hits_<species>.paf`) |
| `insilico_pcr.py` | In silico PCR with the panel primers in each genome → `amplicons_ref.fa` (60 references; MSA sequence used where no primer site was found, see `insilico_pcr.tsv`) |
| `run_sample.sh` | Per sample: 5,000 reads (seqtk, seed 100), minimap2 to all 60 references (all secondary alignments), Kraken2 per-read classification |
| `summarise.py` | Alignment-identity and Kraken2 read-vote species calls; per-locus and leave-one-locus-out accuracy (`summary.txt`, `per_locus_calls.tsv`) |
| `identity_profile.py`, `foreign_reads.py`, `coustani_detail.py` | Read identity to each species' reference; reads closely matching another panel species (Table S4) |
| `bracken900/` | Species calls with Bracken k-mer distributions rebuilt for 900 bp reads (`bracken-build -k 35 -l 900`) |

### 10. Quality control and run checks: `qc_checks/`

Specimens were sequenced in 8 runs on 5 R10.4.1 flow cells, some of which were washed and reused.
These scripts trace every file to its run, barcode and flow cell, and apply one quality-control rule to all species.

**QC rule.** A specimen is kept if it has at least 1,000 reads and is neither
(a) **background**: reads ≤ 1.5× the largest unassigned-barcode read count of its run *and* species
composition within L1 < 15 of the run's pooled composition (`background_check.py`), nor
(b) **carry-over**: reads ≤ the count expected from the previous library on the same flow cell and barcode
(1% of the previous run plus 0.1% of older runs; `carryover.py`, `carryover_flags.py`).
The excluded files and the reason for each are listed in `excluded_qc.tsv` (plus two *An. funestus* files removed as carry-over; see `carryover_flags.tsv`).

| Script | Purpose | Output |
|---|---|---|
| `samples_map.py`, `run_summary.sh` | Specimen → run, barcode, flow cell, read count | `study_samples.tsv`, `runs.tsv`, `verified_runs.tsv` |
| `verify_anoph_speciation.py`, `verify_slca003.py` | Read-ID trace of files back to raw MinKNOW output; on-target and IR-panel checks for SLCA003 | `anoph_speciation_trace.tsv`, `slca003_verification.tsv`, `on_target_fraction.tsv`, `ir_panel_reads.tsv` |
| `subsample.py`, `classify.sh`, `composition.py` | Kraken2 composition of every barcode in every run, including unassigned barcodes | `composition.tsv`, `raw_barcode_counts.tsv`, `slca020_021/` |
| `carryover.py`, `carryover_flags.py` | Same-barcode carry-over between consecutive libraries on a reused flow cell | `carryover_tests.tsv`, `carryover_per_sample.tsv`, `carryover_flags.tsv` |
| `background_check.py` | Cross-barcode background test | `background_runs.tsv`, `background_flags.tsv` |
| `per_run_accuracy.py`, `species_checks.py` | Accuracy by run and per-species checks | `per_run_accuracy.tsv`, `<species>_checks.tsv` |
| `coustani_identity.py`, `novelty_test.py` | Read identity to the 60 amplicon references; read-identity (<95%) "unresolved" rule | `coustani_identity.tsv`, `funestus_identity.tsv`, `novelty_per_sample.tsv` |
| `rerun_after_qc.sh`, `rerun_geno_tree.sh` | Reruns classification, accuracy, depth, genotyping and the tree on the 349-sample set | |
| `build_additional_file2.py` | Additional file 2 (per-specimen metadata and QC workbook) | |

`verified_set_MANIFEST.tsv` lists each analysed file with its raw source run and read-count check.
`manuscript_numbers_349.txt` collects the numbers quoted in the manuscript.

### 11. Study-site map: `figure_scripts/fig1_map/` (Fig. 1)

County boundaries are geoBoundaries gbOpen KEN ADM1 (public domain). Points are the mean GPS position of each
collection village (`collection_sites_clean.tsv`); they are not linked to individual sequenced specimens.

### 12. SRA submission: `sra_submission/`

`build_sra_metadata.py` writes the BioSample (Invertebrate.1.0) and SRA metadata sheets for all 374 sequenced
specimens from `files_manifest.tsv`; `stage_upload_files.sh` gathers the FASTQs as `<sample_id>.fastq.gz` and checks them against `md5_source_files.txt`.

## Resources

- `nomadsSID.amplicons.bed`: target amplicon coordinates on AgamP4.
- `resources/primers.tsv`: primer sequences for the 10 species-ID loci (Table 2).

## License

MIT. See [LICENSE](LICENSE).
