import pandas as pd
import glob

fs = glob.glob("*.json")

l = []
for i, f in enumerate(fs):
    if i % 50 == 0:
        print(f"Processed {i} files...")

    df = pd.read_json(f)

    l.append(df)

df = pd.concat(l)

df.to_json("consolidated.json")
