This is a Snakemake workflow for running [SINGER](https://github.com/popgenmethods/SINGER) (MCMC sampling of ancestral recombination graphs) in parallel (e.g. across chunks of sequence). It largely cannibalizes code from python scripts in the main SINGER repo, with some modification to adjust the parameters for each chunk for missing sequence and recombination rate heterogeneity. Some diagnostic plots are produced at the end, that compare summary statistics to their expectations given the ARG topology.

### Dependencies

Using `git` and `mamba` and `pip`:

```bash
git clone https://github.com/nspope/singer-snakemake my-singer-run && cd my-singer-run
mamba install -c bioconda snakemake
python3 -m pip install -r requirements.txt
snakemake --cores=20 --configfile=configs/example_config.yaml
```

### Inputs

The input files for each chromosome are:

  - __chromosome_name.vcf__ VCF that can be used as SINGER input (diploid, phased, not compressed)
  - __chromosome_name.mask.bed__ (optional) bed file containing inaccessible intervals
  - __chromosome_name.hapmap__ (optional) recombination map in the format described in the documentation for `msprime.RateMap.read_hapmap` (see [here](https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap))

see `example/*`.

### Config

A template for the configuration file is in `configs/example_config.yaml`:

```yaml
# --- example_config.yaml ---
input-dir: "example" # directory with input files per chromosome, that are "chrom.vcf" "chrom.hapmap" "chrom.mask.bed"
max-resumes: 1000 # maximum number of times to try to resume MCMC on error per iteration
chunk-size: 1e6 # target size in base pairs for each singer run
mutation-rate: 1e-8 # per base per generation mutation rate
recombination-rate: 1e-8 # per base per generation recombination rate, ignored if hapmap is present
max-missing: 0.975 # ignore chunks with more than this proportion of missing bases
mcmc-samples: 100 # number of MCMC samples (each sample is a tree sequence)
mcmc-thin: 10 # thinning interval between MCMC samples
mcmc-burnin: 0.2 # proportion of initial samples discarded when computing plots of statistics
polarised: True # are variants polarised so that the reference state is ancestral
random-seed: 1 # random seed
```

### Outputs

The output files for each chromosome will be generated in `results/chromosome_name`:

  - __chromosome_name.adjusted_mu.p__ : `msprime.RateMap` containing adjusted mutation rates (`proportion_accessible_bases * mutation_rate`) in each chunk
  - __chromosome_name.replicate.trees__ : a tree sequence replicate generated by SINGER
  - __chromosome_name.replicate.stats.p__ : "fitted values" for summary statistics (e.g. branch-mode statistics calculated with tskit)
  - __chromosome_name.stats.p__ : "observed values" for summary statistics (e.g. calculated from with `scikit-allel`)
  - __chunks/*__ the raw SINGER output and logs
  - __diversity-trace.png__, __tajima-d-trace.png__ : MCMC trace for fitted nucleotide diversity and Tajima's D
  - __diversity-scatter.png__, __tajima-d-scatter.png__ : Observed vs fitted nucleotide diversity and Tajima's D, across chunks
  - __diversity-skyline.png__, __tajima-d-skyline.png__ : Observed and fitted nucleotide diversity and Tajima's D, across genome position
  - __folded-afs.png__, __unfolded-afs.png__ : Observed vs fitted site frequency spectra
  - __site-density.png__ : Sanity check showing proportion of missing data, proportion variant bases (out of accessible bases), recombination rate across genome position.
