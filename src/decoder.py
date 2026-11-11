import numpy as np
from itertools import combinations, product
from math import isqrt, ceil, log2, sqrt, comb

class LUTDecoder:

    def __init__(self, H, q, max_candidates=200_000):

        self.q = int(q)

        if self.q < 2 or any(
            self.q % d == 0
            for d in range(2, isqrt(self.q) + 1)
        ):
            raise ValueError("q must be prime")

        self.H = np.asarray(H, dtype=int) % self.q

        if self.H.ndim != 2:
            raise ValueError("H must be a 2D matrix")

        self.M, self.n = self.H.shape

        if self.q ** self.M > max_candidates:
            raise ValueError(
                "LUT requires too many error candidates"
            )

        # Syndrome -> minimum-weight error
        self.table = {}

        # Enumerate all error vectors
        for e in product(
            range(self.q), repeat=self.M
        ):

            s = self.syndrome(e)

            previous = self.table.get(s)

            if previous is None:
                self.table[s] = e

            elif (
                self.hamming_weight(e), e
            ) < (
                self.hamming_weight(previous), previous
            ):
                self.table[s] = e

    @staticmethod
    def hamming_weight(e):
        return sum(v != 0 for v in e)

    def syndrome(self, e):

        e = np.asarray(e, dtype=int)

        if e.shape != (self.M,):
            raise ValueError("Invalid error shape")

        return tuple(
            int(v) for v in (self.H.T @ e) % self.q
        )

    def decode(self, syndrome):
        """
        O(1) dictionary lookup after preprocessing.
        """

        s = tuple(int(v) for v in syndrome)

        if len(s) != self.n or any(
            v < 0 or v >= self.q for v in s
        ):
            raise ValueError("Invalid syndrome")

        return self.table.get(s)


class BWSDDecoder:

    def __init__(self, H, q, degree=1):

        self.q = int(q)

        if self.q < 2 or any(
            self.q % d == 0
            for d in range(2, isqrt(self.q) + 1)
        ):
            raise ValueError("q must be prime")

        self.H = np.asarray(H, dtype=int) % self.q

        if self.H.ndim != 2:
            raise ValueError("H must be a 2D matrix")

        self.M, self.n = self.H.shape
        self.degree = int(degree)

        if not 0 <= self.degree <= self.M:
            raise ValueError("Invalid degree")

    def syndrome(self, e):

        e = np.asarray(e, dtype=int)

        if e.shape != (self.M,):
            raise ValueError("Invalid error shape")

        return tuple(
            int(v) for v in (self.H.T @ e) % self.q
        )

    def error_patterns(self):
        """
        Generate errors with Hamming weight <= degree.
        """

        for k in range(self.degree + 1):

            for positions in combinations(
                range(self.M), k
            ):

                for values in product(
                    range(1, self.q), repeat=k
                ):

                    e = [0] * self.M

                    for i, value in zip(
                        positions, values
                    ):
                        e[i] = value

                    yield tuple(e)

    def decode(self, syndrome):
        """
        Search errors in increasing Hamming weight.
        Return the first error matching the syndrome.
        """

        s = tuple(int(v) for v in syndrome)

        if len(s) != self.n or any(
            v < 0 or v >= self.q for v in s
        ):
            raise ValueError("Invalid syndrome")

        for e in self.error_patterns():

            if self.syndrome(e) == s:
                return e

        return None

class CollisionAwareAmplitudeAggregation: 
    def __init__(self, B, q, degree=1, atol=1e-12):
        self.B = np.asarray(B, dtype=int)
        self.q = int(q)
        self.degree = int(degree)
        self.atol = atol
        self.m, self.n = self.B.shape

        if not 0 <= self.degree <= self.m: 
            raise ValueError("Invalid polynomial degree")
        self.B %= self.q 
        # Binary encoding of F_q  
        self.bits_per_symbol = ceil(log2(self.q))
        self.num_qubits = self.n * self.bits_per_symbol

        # Aggregation results 
        self.raw_amplitudes = {}
        self.normalized_amplitudes = {}
        self.syndrome_errors = {}

        self.num_errors = 0 
        self.normalization = None

    def count_errors(self): 
        '''
        Count errors of Hamming weight <= degree
        '''
        return sum(comb(self.m, k)*(self.q-1)**k for k in range(self.degree + 1))

    def enumerate_errors(self): 
        '''
        Generate all errors with weight <= degree
        '''
        for weight in range(self.degree+1): 
            for positions in combinations(range(self.m), weight): 
                for values in product(range(1, self.q), repeat=weight): 
                    error = np.zeros(self.m, dtype=int)
                    for pos, value in zip(positions, values): 
                        error[pos]=value
                    yield error

    def compute_syndrome(self, error): 
        '''
        Calculate s = B^T e mod q 
        '''
        error = np.asarray(error, dtype=int)
        if error.shape != (self.m, ): 
            raise ValueError("Invalid error dimension")
        if np.any(error<0) or np.any(error >= self.q): 
            raise ValueError("Error value outside F_q")
        return tuple((self.B.T@error%self.q).tolist())

    def aggregate(self, error_amplitudes): 
        pass
