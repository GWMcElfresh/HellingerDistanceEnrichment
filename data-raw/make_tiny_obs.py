#!/usr/bin/env python3
"""Regenerate inst/extdata/tiny_obs.h5ad for the extract-metadata vignette.

Call Python anndata directly. Do not use the CRAN R wrapper, which still
passes dtype= (removed from Python anndata >= 0.11).
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "inst" / "extdata" / "tiny_obs.h5ad"

N_CELLS = 60
N_GENES = 10
SUBJECTS = np.repeat([f"S{i}" for i in range(1, 5)], 15)
GROUPS = np.repeat(["A", "A", "B", "B"], 15)
CATEGORIES = np.array([f"Type{i % 3}" for i in range(N_CELLS)])

obs = pd.DataFrame(
    {"subjectId": SUBJECTS, "category": CATEGORIES, "group": GROUPS},
    index=[f"cell{i}" for i in range(N_CELLS)],
)
var = pd.DataFrame(index=[f"gene{i}" for i in range(N_GENES)])
rng = np.random.default_rng(1)
X = rng.normal(size=(N_CELLS, N_GENES)).astype(np.float32)

adata = ad.AnnData(X=X, obs=obs, var=var)
OUT.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUT, compression=None)
print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
