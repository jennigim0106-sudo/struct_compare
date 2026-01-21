# analysis/schrodinger/parse.py


import pandas as pd
from pathlib import Path


def parse_glide_csv(csv_path: str) -> pd.DataFrame:
    csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(f"Glide CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    return df