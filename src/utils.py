import os 
import numpy as np
import pandas as pd 
from itertools import combinations, product

def read_csv(path:str): 
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Cannot found file {path}")
    df = pd.read_csv(path)
    return df.values, df.columns

def check_collisions(B, q, degree):
    m, n = B.shape
    seen = {}
    total = 0

    for k in range(degree + 1):
        for support in combinations(range(m), k):
            for values in product(range(1, q), repeat=k):
                syndrome = np.zeros(n, dtype=int)

                for i, v in zip(support, values):
                    syndrome += v * B[i]

                s = tuple(syndrome % q)
                error = (support, values)
                total += 1

                if s in seen:
                    print("Collision detected!")
                    print("Syndrome:", s)
                    print("First error:", seen[s])
                    print("Second error:", error)
                    return

                seen[s] = error

    print("Candidate errors:", total)
    print("Unique syndromes:", len(seen))
    print("No collisions detected.")
