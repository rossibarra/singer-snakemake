"""
Test data transformations used in validation workflow.
"""

import numpy as np
import msprime
import pytest
from workflow.validation.utils import repolarise_tree_sequence
from workflow.validation.utils import simulate_mispolarisation


def test_repolarise_tree_sequence():
    ts = msprime.sim_ancestry(
        samples=4, 
        sequence_length=100, 
        population_size=1, 
        recombination_rate=0, 
        random_seed=1,
    )
    ts = msprime.sim_mutations(ts, rate=0.1, random_seed=2)
    biallelic = np.bincount(ts.mutations_site, minlength=ts.num_sites) == 1
    assert np.any(~biallelic)
    mispolarise = np.full_like(biallelic, False)
    mispolarise[~biallelic] = True
    mispolarise[:ts.num_sites // 2] = True
    geno_polar = repolarise_tree_sequence(ts, mispolarise).genotype_matrix()
    geno = ts.genotype_matrix()
    # biallelic sites that are flagged are flipped
    flipped = np.logical_and(mispolarise, biallelic)
    assert np.any(flipped)
    assert np.all(np.logical_xor(geno_polar[flipped], geno[flipped]))
    # all other sites are unchangged
    assert np.all(~np.logical_xor(geno_polar[~flipped], geno[~flipped])) 


def test_major_allele_repolarisation():
    ts = msprime.sim_ancestry(
        samples=10, 
        sequence_length=1e5,
        population_size=1e4, 
        recombination_rate=1e-8, 
        random_seed=1,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=2)
    biallelic = np.bincount(ts.mutations_site, minlength=ts.num_sites) == 1
    ts = ts.delete_sites(~biallelic)
    ts = repolarise_tree_sequence(ts, simulate_mispolarisation(ts, "maf"))
    geno = ts.genotype_matrix()
    assert np.all(geno.sum(axis=-1) / ts.num_samples <= 0.5)

