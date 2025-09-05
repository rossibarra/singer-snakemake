"""
Run POLEGON on a chunk of ARG produced by SINGER.

Part of https://github.com/nspope/singer-snakemake.
"""

import os
import shutil
import yaml
import subprocess
import numpy as np
import msprime
import tszip
import tskit
from datetime import datetime

# --- lib --- #

def tag(): 
    return f"[singer-snakemake::{snakemake.rule}::{str(datetime.now())}]"


# --- implm --- #

use_polegon = True
if use_polegon:
    params = yaml.safe_load(open(snakemake.input.params))
    params = params.pop("polegon")
    seed = params.pop("seed") + int(snakemake.wildcards.rep)

    # This will match what polegon_master does. 
    # Why the factor of 2 relative to SINGER? 
    # Doesn't it cancel during rescaling, regardless?
    # params["Ne"] *= 2  

    # POLEGON expects inputs named slightly differently than SINGER output
    prefix = snakemake.input.muts.replace("_muts_", "_").removesuffix(".txt")
    shutil.copy(snakemake.input.muts, f"{prefix}_muts.txt")
    shutil.copy(snakemake.input.nodes, f"{prefix}_nodes.txt")
    #shutil.copy(snakemake.input.branches, f"{prefix}_branches.txt")

    # Adjust branch spans to reflect missing data, and absorb mutation rate into span.
    # This is necessary because POLEGON doesn't use the "mutation-adjusted spans"
    # from the mutation map during rescaling. Manually editing the branch spans
    # fixes this.
    with open(snakemake.log.out, "w") as log:
       log.write(
           "f{tag()} Adjusting branch spans input ({prefix}_branches.txt) "
           "to reflect mutational rate map and setting mutation rate to unity\n"
       )
    mutation_map = params.pop("mutation_map")
    adjusted_mu = np.loadtxt(mutation_map, ndmin=2)
    assert np.all(adjusted_mu[:-1, 1] == adjusted_mu[1:, 0])
    adjusted_mu = msprime.RateMap(
       position=np.append(adjusted_mu[0, 0], adjusted_mu[:, 1]),
       rate=adjusted_mu[:, 2],
    )
    branches = np.loadtxt(snakemake.input.branches)
    for i in range(2):
       branches[:, i] = adjusted_mu.get_cumulative_mass(branches[:, i])
    np.savetxt(f"{prefix}_branches.txt", branches)
    params["m"] = 1.0

    invocation = [
        f"{snakemake.params.polegon_binary}",
        "-input", f"{prefix}",
        "-seed", f"{seed}",
    ]
    for arg, val in params.items():
        invocation += f"-{arg} {val}".split()
    
    with open(snakemake.log.out, "a") as out, open(snakemake.log.err, "w") as err:
        print(f"{tag()}", " ".join(invocation), file=out, flush=True)
        process = subprocess.run(invocation, check=False, stdout=out, stderr=err)
        print(f"{tag()} POLEGON run ended ({process.returncode})", file=out, flush=True)

    for suffix in ["muts", "branches", "nodes"]:
        os.remove(f"{prefix}_{suffix}.txt")

    assert process.returncode == 0, f"POLEGON terminated with error ({process.returncode})"
    os.rename(f"{prefix}_new_nodes.txt", snakemake.output.nodes)
else:
    shutil.copy(snakemake.input.nodes, snakemake.output.nodes)
