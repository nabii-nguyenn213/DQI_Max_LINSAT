from math import ceil, isqrt, log2
import numpy as np 
import pennylane as qml
from itertools import combinations, product 
from src.decoder import LUTDecoder, BWSDDecoder

class DQI: 
    def __init__(self, B, allowed, weights, q, 
                 degree=1, seed=42, coefficients=None, 
                 max_states=100_000):
        self.q = int(q) # Finite field 
        if self.q < 2 or any(self.q % k == 0 for k in range(2, isqrt(self.q) + 1)):
            raise ValueError("q must be prime")
        self.B = np.asarray(B, dtype=int) % self.q # Constraint matrix
        self.allowed = np.asarray(allowed, dtype=bool) # Acceptance sets 
        self.weights = np.asarray(weights, dtype=float) # Constraint weights

        self.m, self.n = self.B.shape # m: number of constraints, n: number of lecturers
        if self.allowed.shape != (self.m, self.q): 
            raise ValueError("Allowed must have shape (m, q)")

        self.degree = int(degree) # Polynomial degree 
        self.rng = np.random.default_rng(seed)
        self.max_states = max_states

        r = self.allowed.sum(axis=1) # Size of acceptance sets 
        self.active = np.flatnonzero(
                                    (r>0)& # a constraint must accept at least one value 
                                    (r<self.q)& # a constraint must not accept every field value
                                    (self.weights>0)
                                    )
        # For the tiny example, all 6 constaints are satisfied -> self.active = [0, 1, 2, 3, 4, 5]
        self.H = self.B[self.active]
        self.M = len(self.active)
        
        # NOTE : Quantum register
        self.bits = (self.q-1).bit_length() # Number of qubits per field element
        
        self.error_qubits = self.M * self.bits
        self.syndrome_qubits = self.n * self.bits

        self.num_qubits = self.error_qubits + self.syndrome_qubits

        self.error_wires = list(range(self.error_qubits))
        self.syndrome_wires = list(range(self.error_qubits, self.num_qubits))

        self.dev = qml.device("default.qubit", wires=self.num_qubits) # Pennylane simulator

        # Calculate Fourier coefficients 
        self.fourier = self.fourier_coefficients()
        # Construct normalized error layers 
        self.layers = self.build_error_layer()
        # Initialize normalized polynomial coefficients.
        # They will be optimized after selecting a decoder.
        self.coefficients = np.ones(
            self.degree + 1,
            dtype=complex
        )
        self.coefficients /= np.linalg.norm(self.coefficients)
        self.error_state = self.build_error_state()
        self.last_success_probability = None

    def objective(self, y): 
        '''
        Calculate Max-LINSAT objective S(y)
        '''
        y = np.asarray(y, dtype=int)
        values = (self.B @ y) % self.q
        satisfied = self.allowed[np.arange(self.m), values]
        return float(np.dot(self.weights, satisfied))

    def info(self): 
        return {
                'field': self.q, 
                'variables': self.n, 
                'constraints': self.m, 
                'active_constraints': self.M, 
                'degree': self.degree, 
                'error_qubits': self.error_qubits, 
                'syndrome_qubits': self.syndrome_qubits, 
                'total_qubits': self.num_qubits
                }

    def fourier_coefficients(self): 
        '''
        Compute Finite-field Fourier coefficients
        for every active Max-LINSAT constraint.

        Returns:
            Fourier matrix with shape (M, q)

        Each row corresponds to an active constraint.
        Each column corresponds to frequency k.
        '''
        fourier = np.zeros((self.M, self.q), dtype=complex) # Fourier coefficient matrix
        z = np.arange(self.q) # All elements of F_q 

        for t, i in enumerate(self.active): 
            # Acceptance set as an indicator function
            f = self.allowed[i].astype(float)
            # Number of accepted values
            r = int(np.sum(f))
            # Center and normalize
            g = (self.q * f - r)/np.sqrt(r * (self.q - r))
            # Calculate Fourier coefficients
            for k in range(self.q):
                fourier[t, k] = np.sum(g * np.exp(2j * np.pi * k * z / self.q)) / self.q
        return fourier

    def build_error_layer(self): 
        '''
        Construct normalized Fourier error layers 
        for degrees k = 0, ..., l 
        
        Returns: 
            layers[k]={error_vector: complex_amplitude}
        '''
        layers = []
        for k in range(self.degree+1): 
            layer = {}
            # Choose positions of nonzero errors 
            for positions in combinations(range(self.M), k): 
                # Assign values 1, ..., q-1
                for values in product(range(1, self.q), repeat=k):
                    e = [0]*self.M 
                    amplitude = 1.0+0j 
                    for i, v in zip(positions, values):
                        e[i] = v
                        # Fourier coefficient and
                        # weight of active constraint i
                        amplitude *= self.weights[self.active[i]] * self.fourier[i, v]
                    layer[tuple(e)] = amplitude
            # Normalize the entire degree-k layer
            norm = np.sqrt(sum(abs(a)**2 for a in layer.values()))
            if norm < 1e-12:
                raise ValueError(f"Degree {k} has zero amplitude")
            layer = {e: a / norm for e, a in layer.items()}
            layers.append(layer)
        return layers

    def polynomial_coefficients(self, decoder):
        """
        Optimize degree-l polynomial coefficients for
        the postselected output of a selected decoder.

        Suitable for small instances using exact
        classical enumeration.
        """
        assignments = list(
            product(range(self.q), repeat=self.n)
        )
        N = len(assignments)
        if N > 100_000:
            raise ValueError(
                "Instance too large for exact optimization"
            )
        # Each column corresponds to one polynomial degree
        basis = np.zeros(
            (N, self.degree + 1),
            dtype=complex
        )
        for k, layer in enumerate(self.layers):
            state = np.zeros(
                (self.q,) * self.n,
                dtype=complex
            )
            for e, amplitude in layer.items():

                s = self.syndrome(e)
                decoded = decoder.decode(s)
                # Only successfully decoded branches
                # survive postselection.
                if decoded is not None and tuple(decoded) == e:
                    state[s] += amplitude
            basis[:, k] = np.fft.fftn(
                state,
                norm="ortho"
            ).ravel()
        # Objective values for all assignments
        scores = np.array([
            self.objective(y)
            for y in assignments
        ])
        # Gram matrix
        G = basis.conj().T @ basis
        # Projected objective matrix
        K = basis.conj().T @ (
            scores[:, None] * basis
        )
        G = (G + G.conj().T) / 2
        K = (K + K.conj().T) / 2
        # Generalized eigenvalue optimization.
        # Handle zero-norm degree layers safely.
        eigenvalues_G, eigenvectors_G = np.linalg.eigh(G)
        valid = eigenvalues_G > 1e-10
        if not np.any(valid):
            raise ValueError(
                "No error states survive decoding"
            )
        W = (
            eigenvectors_G[:, valid]
            / np.sqrt(eigenvalues_G[valid])
        )
        K_reduced = W.conj().T @ K @ W
        K_reduced = (
            K_reduced + K_reduced.conj().T
        ) / 2
        eigenvalues, eigenvectors = np.linalg.eigh(
            K_reduced
        )
        coefficients = W @ eigenvectors[:, -1]
        # Normalize the original error superposition
        coefficients /= np.linalg.norm(coefficients)
        # Remove global phase
        first = np.flatnonzero(
            np.abs(coefficients) > 1e-12
        )[0]
        coefficients *= np.exp(
            -1j * np.angle(coefficients[first])
        )
        self.K = K
        self.G = G
        self.projected_score = float(eigenvalues[-1])
        return np.real_if_close(coefficients)

    def build_error_state(self): 
        """
        Construct the normalized error-state vector.

        Returns:
            state: Complex vector of size
                2 ** self.error_qubits
        """
        # Number of computational basis states
        dimension = 2 ** self.error_qubits
        # Initialize all amplitudes to zero
        state = np.zeros(
            dimension,
            dtype=complex
        )
        # Iterate over polynomial degrees
        for k, layer in enumerate(self.layers):
            # Polynomial coefficient for degree k
            w_k = self.coefficients[k]
            # Iterate over the error vectors
            for e, amplitude in layer.items():
                # Convert error vector into binary index
                index = self.encode_error(e)
                # Store the final amplitude
                state[index] = w_k * amplitude
        # Verify normalization
        norm = np.linalg.norm(state)
        if not np.isclose(norm, 1.0):
            raise ValueError(
                f"Error state is not normalized: {norm}"
            )
        return state

    def prepare_error_state(self): 
        """
        Prepare the quantum error register using PennyLane.

        Must be called inside a PennyLane QNode.
        """

        qml.StatePrep(
            self.error_state,
            wires=self.error_wires
        )

    def syndrome(self, e): 
        '''
        Compute the syndrome H.T @ e mod q
        '''
        e = np.asarray(e, dtype=int)
        if e.shape != (self.M,):
            raise ValueError("Invalid error vector shape")
        s = (self.H.T @ e) % self.q
        return tuple(int(v) for v in s)

    def encode_error(self, e): 
        '''
        Convert an F_q error vector into a binary 
        computational basis index 
        '''
        index = 0 
        for value in e: 
            index = index*(2**self.bits)+int(value)
        return index

    def create_decoder(self, decoder_type="bwsd"):
        if decoder_type == "lut":
            return LUTDecoder(
                H=self.H,
                q=self.q
            )
        elif decoder_type == "bwsd":
            return BWSDDecoder(
                H=self.H,
                q=self.q,
                degree=self.degree
            )
        else:
            raise ValueError(
                "decoder_type must be 'lut' or 'bwsd'"
            )

    def verify_decoder(self, decoder):
        """
        Validate decoder corrections.

        Syndrome collisions are allowed.
        Returns syndrome -> selected error.
        """
        if decoder.q != self.q or not np.array_equal(
            decoder.H, self.H
        ):
            raise ValueError(
                "Decoder configuration does not match DQI"
            )
        corrections = {}
        for layer in self.layers:
            for e in layer:
                s = self.syndrome(e)
                if s in corrections:
                    continue
                decoded = decoder.decode(s)
                if decoded is None:
                    continue
                decoded = tuple(int(v) for v in decoded)
                if len(decoded) != self.M:
                    raise ValueError(
                        "Invalid decoded error length"
                    )
                if self.syndrome(decoded) != s:
                    raise ValueError(
                        f"Incorrect decoder output for {s}"
                    )
                corrections[s] = decoded
        return corrections

    def apply_decoder(self, decoder):
        """
        Reversibly remove the error information:

        |e>|s> -> |e - D(s)>|s>

        Must be called inside a PennyLane QNode.
        """
        # Validate decoder on all supported errors
        corrections = self.verify_decoder(decoder)
        # Group qubits into field-element registers
        E = [
            self.error_wires[
                i * self.bits : (i + 1) * self.bits
            ]
            for i in range(self.M)
        ]
        # Precompute modular subtraction gates
        subtraction_gates = {
            value: self.subtract_unitary(value)
            for value in range(1, self.q)
        }
        # Apply syndrome-controlled corrections
        for syndrome, error in corrections.items():
            # Encode syndrome into binary control bits
            control_bits = []
            for value in syndrome:
                binary = format(
                    value,
                    f"0{self.bits}b"
                )
                control_bits.extend(
                    int(bit) for bit in binary
                )
            # Subtract the decoded error
            for i, value in enumerate(error):
                if value == 0:
                    continue
                qml.ControlledQubitUnitary(
                    subtraction_gates[value],
                    wires=list(self.syndrome_wires) + E[i],
                    control_values=control_bits
                )

    def decode(self): 
        '''
        Connect the decoder 
        '''

    def subtract_unitary(self, value):
        """
        |x> -> |x - value mod q>
        """
        d = 2 ** self.bits
        U = np.zeros((d, d), dtype=complex)
        for x in range(d):
            if x < self.q:
                new_x = (x - value) % self.q
            else:
                new_x = x
            U[new_x, x] = 1.0
        return U

    def modular_add_unitary(self, coefficient):
        """
        Construct a reversible modular addition gate.

        |e>|s> -> |e>|s + coefficient*e mod q>
        """
        q = self.q
        d = 2 ** self.bits
        # Matrix dimension: (d*d) x (d*d)
        U = np.zeros(
            (d * d, d * d),
            dtype=complex
        )
        coefficient = int(coefficient) % q
        for e in range(d):
            for s in range(d):
                # Valid finite-field states
                if e < q and s < q:
                    new_s = (s + coefficient * e) % q
                # Leave unused binary states unchanged
                else:
                    new_s = s
                # Input computational basis index
                input_index = e * d + s
                # Output computational basis index
                output_index = e * d + new_s
                U[output_index, input_index] = 1.0
        return U

    def compute_syndrome(self):
        """
        Quantum syndrome computation:

        |e>|s> -> |e>|s + H.T @ e mod q>

        When s starts at zero:
        |e>|0> -> |e>|H.T @ e>
        
        Must be called inside a PennyLane QNode.
        """
        # Divide error qubits into field-element registers
        E = [
            self.error_wires[
                i * self.bits : (i + 1) * self.bits
            ]
            for i in range(self.M)
        ]
        # Divide syndrome qubits into field-element registers
        S = [
            self.syndrome_wires[
                j * self.bits : (j + 1) * self.bits
            ]
            for j in range(self.n)
        ]
        # Cache the modular addition gates
        addition_gates = {}
        for c in np.unique(self.H):
            c = int(c)
            if c != 0:
                addition_gates[c] = (
                    self.modular_add_unitary(c)
                )
        # Apply syndrome computation
        for i in range(self.M):
            for j in range(self.n):
                c = int(self.H[i, j])
                if c == 0:
                    continue
                qml.QubitUnitary(
                    addition_gates[c],
                    wires=E[i] + S[j]
                )

    def syndrome_circuit(self):
        """
        Prepare the DQI error superposition
        and compute the syndrome quantum mechanically.
        """
        @qml.qnode(self.dev)
        def circuit():
            # Step 4: Error superposition
            self.prepare_error_state()
            # Step 5: Syndrome computation
            self.compute_syndrome()
            # Syndrome measurement probabilities
            return qml.probs(
                wires=self.syndrome_wires
            )
        return circuit()

    def inverse_qft_unitary(self):
        """
        Construct the inverse Fourier transform over F_q.

        Unused binary states are left unchanged.
        """
        q = self.q
        d = 2 ** self.bits
        # Identity matrix for the full qubit register
        U = np.eye(d, dtype=complex)
        # Finite-field values
        x = np.arange(q)
        # Construct inverse QFT over F_q
        F_inverse = np.exp(
            -2j * np.pi * np.outer(x, x) / q
        ) / np.sqrt(q)
        # Embed into the qubit register
        U[:q, :q] = F_inverse
        return U

    def apply_inverse_qft(self):
        """
        Apply inverse QFT to the syndrome register.

        Must be called inside a PennyLane QNode.
        """
        # Construct inverse Fourier unitary
        U = self.inverse_qft_unitary()
        # Apply it to each syndrome field element
        for j in range(self.n):
            wires = self.syndrome_wires[
                j * self.bits : (j + 1) * self.bits
            ]
            qml.QubitUnitary(
                U,
                wires=wires
            )

    def output_probabilities(self, decoder_type="bwsd"):
        """
        Execute PennyLane DQI with postselected decoding.

        Returns assignment probabilities conditioned on
        the error register being all zero.
        """
        # 1. Select decoder
        decoder = self.create_decoder(decoder_type)
        # 2. Optimize polynomial for this decoder
        self.coefficients = self.polynomial_coefficients(
            decoder
        )
        # 3. Rebuild the error state
        self.error_state = self.build_error_state()
        # 4. Validate corrections
        self.verify_decoder(decoder)
        # 5. Construct PennyLane QNode
        @qml.qnode(self.dev)
        def circuit():
            # Prepare bounded-weight error superposition
            self.prepare_error_state()
            # Compute H.T @ e
            self.compute_syndrome()
            # Reversible decoder correction
            self.apply_decoder(decoder)
            # Inverse QFT on syndrome only
            self.apply_inverse_qft()
            # Measure both registers
            return qml.probs(
                wires=(
                    self.error_wires
                    + self.syndrome_wires
                )
            )
        raw = np.asarray(circuit(), dtype=float)
        # 6. Select error = |0...0> branch
        syndrome_dimension = 2 ** len(
            self.syndrome_wires
        )
        # Error wires come first in the measurement order.
        successful = raw[:syndrome_dimension]
        success_probability = float(
            np.sum(successful)
        )
        self.last_success_probability = success_probability
        if success_probability < 1e-12:
            raise ValueError(
                "No successful decoding branches"
            )
        # 7. Conditional syndrome probabilities
        successful /= success_probability
        # 8. Convert binary basis indices to F_q assignments
        d = 2 ** self.bits
        probabilities = {}
        for y in product(
            range(self.q),
            repeat=self.n
        ):
            index = 0
            for value in y:
                index = index * d + value
            probabilities[y] = float(
                successful[index]
            )
        total = sum(probabilities.values())
        if not np.isclose(total, 1.0, atol=1e-8):
            raise ValueError(
                "Invalid postselected distribution"
            )
        return probabilities

    def reference_probabilities(self, decoder):
        """
        Classical reference for postselected DQI.

        Assumes self.coefficients have already been
        optimized for the selected decoder.
        """
        state = np.zeros(
            (self.q,) * self.n,
            dtype=complex
        )
        for k, layer in enumerate(self.layers):
            for e, amplitude in layer.items():
                s = self.syndrome(e)
                decoded = decoder.decode(s)
                if decoded is not None and tuple(decoded) == e:
                    state[s] += (
                        self.coefficients[k] * amplitude
                    )
        output = np.fft.fftn(
            state,
            norm="ortho"
        )
        probabilities = np.abs(output) ** 2
        success = float(np.sum(probabilities))
        if success < 1e-12:
            raise ValueError("Postselection probability is zero")
        probabilities /= success
        return {
            y: float(probabilities[y])
            for y in np.ndindex(probabilities.shape)
        }

    def expected_score(self, probabilities):
        """
        Expected Max-LINSAT score under a
        given assignment probability distribution.
        """
        return sum(
            probability * self.objective(y)
            for y, probability in probabilities.items()
        )

    def circuit(self): 
        '''
        Build and run the Pennylane circuit
        '''
