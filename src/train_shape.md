
#!/usr/bin/env python3
"""
# --- *_reps.npz ---
# Key: 'embeddings     ' | Shape: (19667, 960)         | Dtype: float32
# Key: 'ids            ' | Shape: (19667,)             | Dtype: <U10
"""

from pyprojroot import here
import torch as pt
from pathlib import Path
import numpy as np

DIR_REPS = here() / 'data/reps/'
npz_files = list(DIR_REPS.glob('*.npz'))
print(f"Found {len(npz_files)} .npz files.")

