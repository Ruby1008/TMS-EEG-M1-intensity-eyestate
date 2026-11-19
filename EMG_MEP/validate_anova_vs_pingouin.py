"""Cross-check rm_anova_2x3 against pingouin on random and on test data."""
import sys, os
import numpy as np, pandas as pd, pingouin as pg
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mep_reprocess import rm_anova_2x3, CELLS

rng = np.random.default_rng(1)
Y = rng.normal(size=(15, 6)) + np.array([0, .3, .9, .1, .4, 1.2]) + rng.normal(size=(15, 1))
Y[:, 2] += rng.normal(scale=1.5, size=15)
mine = rm_anova_2x3(Y)
long = pd.DataFrame([(s, e, i, Y[s, k]) for s in range(15) for k, (e, i) in enumerate(CELLS)], columns=["s", "eye", "int", "y"])
ref = pg.rm_anova(long, dv="y", within=["eye", "int"], subject="s", correction=True)
print(mine.round(4).to_string())
print(ref.round(4).to_string())
