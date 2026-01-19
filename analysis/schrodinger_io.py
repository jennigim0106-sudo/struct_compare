import os
import zipfile
import pandas as pd

def extract_schrodinger_zip(zip_path, out_dir):
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(out_dir)
    print("=== WALK ROOT ===", out_dir)
    for root, dirs, files in os.walk(out_dir):
        print(root)
        for f in files:
            print("  ", f)

    return out_dir


def parse_glide_csv(extract_dir):
    for root, _, files in os.walk(extract_dir):
        for f in files:
            if f.endswith(".csv") and "glide" in f.lower():
                csv_path = os.path.join(root, f)
                df = pd.read_csv(csv_path)

                return df.rename(columns={
                    "i_i_glide_lignum": "pose_id",
                    "r_i_glide_gscore": "glide_score"
                })[["pose_id", "glide_score"]]

    raise FileNotFoundError("Glide CSV not found in Schrödinger zip")

