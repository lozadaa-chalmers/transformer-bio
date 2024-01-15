import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--save-path")

args = parser.parse_args()

d = {'col1': [1, 2], 'col2': [3, 4]}
df = pd.DataFrame(data=d)

df.to_csv(args.save_path + "/toy_dataframe.csv", index=False)

df_new = pd.read_csv(args.save_path + "/toy_dataframe.csv")

