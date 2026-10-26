import numpy as np
from itertools import combinations, product
from math import isqrt

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
