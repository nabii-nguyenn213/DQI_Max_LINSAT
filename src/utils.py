import os 
import numpy as np
import pandas as pd 

def read_csv(path:str): 
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Cannot found file {path}")
    df = pd.read_csv(path)
    return df.values, df.columns
