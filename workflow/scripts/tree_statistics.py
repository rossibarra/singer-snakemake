"""
Calculate branch statistics from a tree sequence.  This is done via simulation,
so as to account for missing data without excessive windowing.

Part of https://github.com/nspope/singer-snakemake.
"""

import msprime
import tskit
import tszip
import pickle
import numpy as np
import yaml
import itertools
from collections import defaultdict
from datetime import datetime


# --- lib --- #

def tag(): 
    return f"[singer-snakemake::{snakemake.rule}::{str(datetime.now())}]"


# --- implm --- #

windows = pickle.load(open(snakemake.input.windows, "rb"))
adjusted_mu = pickle.load(open(snakemake.input.mut_rate, "rb"))
inaccessible = pickle.load(open(snakemake.input.inaccessible, "rb"))
seed = 1 + int(snakemake.wildcards.rep) * 1000

accessible = msprime.RateMap(
    position=inaccessible.position,
    rate=1 - inaccessible.rate,
)
accessible_bp = np.diff(accessible.get_cumulative_mass(windows.position))

ts = msprime.sim_mutations(
    tszip.decompress(snakemake.input.trees), 
    rate=adjusted_mu, 
    random_seed=seed, 
    keep=False,
)

diversity = ts.diversity(
    mode='site', 
    windows=windows.position, 
    span_normalise=False,
)
diversity[windows.rate == 1.0] /= accessible_bp[windows.rate == 1.0]
diversity[windows.rate == 0.0] = np.nan

tajima_d = ts.Tajimas_D(
    mode='site', 
    windows=windows.position, 
)
tajima_d[windows.rate == 0.0] = np.nan

afs = ts.allele_frequency_spectrum(
    mode='site', 
    span_normalise=False,
    polarised=snakemake.params.polarised,
) / accessible_bp[windows.rate == 1.0].sum()

stats = {
    "diversity": diversity, 
    "tajima_d": tajima_d, 
    "afs": afs, 
}
pickle.dump(stats, open(snakemake.output.stats, "wb"))

# TODO: add proportion multimapped sites, to use as a trace diagnostic

# TODO: populations are now encoded directly, so use ts.samples(population=i)
# stratified summary stats
strata_stats = {}
if snakemake.params.stratify is not None:
    sample_sets = defaultdict(list)
    for ind in ts.individuals():
        strata = ind.metadata[snakemake.params.stratify]
        sample_sets[strata].extend(ind.nodes)
    strata = [n for n in sample_sets.keys()]
    strata_stats["strata"] = strata
    indexes = list(itertools.combinations_with_replacement(range(len(strata)), 2))

    strata_divergence = ts.divergence(
        sample_sets=[s for s in sample_sets.values()], 
        indexes=indexes,
        mode='site',
        windows=windows.position,
        span_normalise=False,
    )
    strata_divergence[windows.rate == 1.0] /= \
        accessible_bp[windows.rate == 1.0, np.newaxis]
    strata_divergence[windows.rate == 0.0] = np.nan
    strata_stats["divergence"] = strata_divergence

    strata_afs = []
    for strata, sample_set in sample_sets.items():
        afs = ts.allele_frequency_spectrum(
            sample_sets=[sample_set], 
            mode='site', 
            polarised=snakemake.params.polarised,
            span_normalise=False,
        ) / accessible_bp[windows.rate == 1.0].sum()
        strata_afs.append(afs)
    strata_stats["afs"] = strata_afs

pickle.dump(strata_stats, open(snakemake.output.strata_stats, "wb"))
