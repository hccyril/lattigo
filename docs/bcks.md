Here is the comprehensive, software-focused implementation development manual based on the provided cryptographic paper "Bootstrapping Bits with CKKS".

# Lattigo-Based Implementation Development Manual: Bootstrapping Bits with CKKS

## 1. Core Document Purpose

This manual is tailored for a software engineer transitioning into cryptography doctoral research, to translate the academic paper’s abstract cryptographic design ("Bootstrapping Bits with CKKS") into actionable, implementation-ready specifications aligned with the `Lattigo` lattice cryptography codebase. It prioritizes software-focused abstraction of mathematical concepts, step-by-step implementation roadmaps, and testable deliverables to bypass steep math learning curves while fully replicating the paper’s functional logic (with reduced parameter dimensions for prototyping efficiency).

---

## 2. Paper Scheme Comprehensive Explanation (Software Engineer-Focused)

### 2.1 Paper Contextual Summary

#### **Background**: 
Fully Homomorphic Encryption (FHE) allows arbitrary computations directly on encrypted data. Traditionally, schemes like DM/CGGI are optimized for binary operations (bits), while CKKS is tailored for real and complex numbers. However, recent advances showed that CKKS can be used for binary operations by treating a bit as a real number close to 0 or 1, providing massive throughput advantages due to its Single-Instruction Multiple-Data (SIMD) batching.


#### **Identified Problem** 
When using standard CKKS for binary data, the "bootstrapping" process (which resets the noise in a ciphertext so computations can continue) is highly inefficient. Standard CKKS requires a massive gap (usually 10 bits) between the base modulus $q_0$ and the scaling factor $\Delta_0$ to accurately approximate the modular reduction function. This consumes an excessive amount of the "modulus budget," leaving less room for actual program logic. Furthermore, standard CKKS sine-based bootstrapping only linearly reduces noise, requiring additional costly "cleaning" functions.

#### **Core Contributions & Innovation Points**:
1. **BinBoot (Binary Bootstrapping)**: 
A novel algorithmic optimization that replaces the standard CKKS sine evaluation with a custom trigonometric function, $f_{BinBoot}$. This allows the scaling factor to be nearly equal to the modulus ($q_0 = 2\Delta_0$), saving over 100 bits of modulus budget while inherently providing quadratic noise reduction.

2. **GateBoot (Combined Bootstrapping and Logic)**: 
A new primitive that executes a binary gate (like NAND or XOR) *simultaneously* with the bootstrapping operation by processing the algebraic sum of two input ciphertexts.

3. **Cross-Scheme Compatibility**: 
The format of `GateBoot` relies on a modulus ratio of $q_0 = 3\Delta_0$, allowing for direct integration with packed DM/CGGI LWE ciphertexts via ring packing techniques like HERMES.


#### **Scheme Full Functional Description (Data Flow)**:

Input plaintext binary array 
$\rightarrow$ Encode as CKKS slots with low scaling factor 
$\rightarrow$ Encrypt 
$\rightarrow$ Perform SIMD binary logic (AND, OR, XOR via real arithmetic) 
$\rightarrow$ When noise threshold is reached, apply `BinBoot` (to clean single ciphertext) or `GateBoot` (to evaluate a gate and clean two ciphertexts) 
$\rightarrow$ Decrypt 
$\rightarrow$ Decode back to binary array.

### 2.2 Algorithm Step-by-Step Breakdown

**1. Cryptographic Parameter Setup Module**

**Input**: Security level $\lambda$, ring dimension $N=2^{14}$ (reduced parameter Param14), base scaling factor $\Delta_0$, default scaling factor $\Delta$.
 
**Output**: CKKS parameter set (`params.Parameters`) where $q_0 = 2\Delta_0$ (for BinBoot) or $q_0 = 3\Delta_0$ (for GateBoot).

**2. BinBoot (Algorithm 2)**
 
**Preconditions**: Input ciphertext $ct$ encrypting $\varphi \in \{0, 1\}^{N/2}$ with accumulated noise.
 
**Step 1**: Execute `CtS` (Coefficient-to-Slots) $\circ$ `ModRaise` $\circ$ `StC` (Slots-to-Coefficients).
 
**Step 2**: Extract the real part homomorphically: $ct^{\prime\prime} = (conj(ct^\prime) + ct^\prime) / 2$.
 
**Step 3**: Homomorphically evaluate the custom trigonometric approximation $Eval_{f_{BinBoot}}(ct^{\prime\prime})$.
 
**Postconditions**: Returned $ct_{out}$ decrypts to the cleaned binary vector.

**3. GateBoot (Algorithm 3)**

**Preconditions**: Two input ciphertexts $ct_1, ct_2$ encrypting binary vectors $\varphi_1, \varphi_2$.
 
**Step 1**: Integer addition: $ct_{add} = ct_1 + ct_2$ (plaintext domain $\{0, 1, 2\}$).
 
**Step 2**: Execute `CtS` $\circ$ `ModRaise` $\circ$ `StC`.
 
**Step 3**: Extract the real part homomorphically.
 
**Step 4**: Homomorphically evaluate the specific gate polynomial approximation $Eval_{f_G}(ct^{\prime\prime})$ (e.g., NAND, AND).
 
**Postconditions**: Returned $ct_{out}$ decrypts to the result of the binary gate applied to $\varphi_1$ and $\varphi_2$.


### 2.3 Full Formula Documentation with Software-Centric Explanations
 
**BinBoot Base Ratio**: $q_0 = 2\Delta_0$ 

**Software Impact**: Configures the Lattigo parameter generation. The base prime $q_0$ must be set to exactly twice the value of the scaling factor, drastically shrinking the standard memory footprint required for LWE/RLWE arrays.

**GateBoot Base Ratio**: $q_0 = 3\Delta_0$ 

**Software Impact**: Configures the Lattigo parameter generation for GateBoot tests. Requires the base prime to be three times the scaling factor to provide equally spaced plaintext points at $0, 1/3, 2/3$.

**BinBoot Evaluation Function**: $f_{BinBoot}(x) = \frac{1-\cos(2\pi x)}{2}$ 

**Software Impact**: This formula represents the mathematical continuous function we must approximate. In the software, this translates to feeding these bounds and target curves into Lattigo's Chebyshev polynomial approximation generator (replacing the default sine function used in CKKS `EvalMod`).

**GateBoot NAND Function**: $f_{NAND}(x) = \frac{2}{3}(1 + \sin(2\pi x + \frac{\pi}{6}))$ 

**Software Impact**: Used when dynamically generating the evaluation polynomial for a SIMD NAND gate. It maps $0 \rightarrow 1$, $1/3 \rightarrow 1$, and $2/3 \rightarrow 0$.


## 3. Lattigo-Based Implementation Specification

### 3.1 Core Implementation Constraints

**Reduced Parameter Target (Param14)**: Ring dimension $N=2^{14}$, sparse secret key Hamming weight $h=32$, dense Hamming weight $h=256$, and maximum key switching depth $dnum=13$. This ensures execution on consumer hardware within the 5-minute threshold.

**Integration Principle**: Leverage standard Lattigo `ckks.Bootstrapper` interfaces but inject custom polynomial evaluation structures during the `EvalMod` phase.

### 3.2 New Entry Point Creation

Create the following structure in the Lattigo fork:

- `/examples/paper-scheme/main.go`: Orchestrates the flow. Initializes `Param14`, generates keys, encrypts two binary arrays, calls the custom `GateBootNAND` function, and asserts the decrypted output matches standard bitwise NAND.
- `/examples/paper-scheme/parameters.go`: Hardcodes the `Param14` modulus chain, specifically enforcing $q_0 = 3\Delta_0$ for GateBoot compatibility.
- `/examples/paper-scheme/utils.go`: Test vector generators (random binary arrays mapped to float64 `[0.0, 1.0]`) and error-logging validators.
- `/examples/paper-scheme/README.md`: Execution instructions (`go run .`) and expected execution times (target: ~1.39s per bootstrap for $N=2^{14}$).

### 3.3 Reusable Lattigo Module Methods

1. `github.com/tuneinsight/lattigo/v5/core/rlwe`: For core KeyGenerator instantiation (`NewKeyGenerator`), Encryptor (`NewEncryptor`), and Decryptor (`NewDecryptor`).
2. `github.com/tuneinsight/lattigo/v5/ckks`: For standard CKKS encoding (`NewEncoder`), decoding, and basic addition operations (`Evaluator.Add`).
3. `github.com/tuneinsight/lattigo/v5/ckks/bootstrapping`: For invoking standard `SlotsToCoeffs` (StC) and `CoeffsToSlots` (CtS) transformations without rewriting the linear transformations.
4. `github.com/tuneinsight/lattigo/v5/utils/bignum`: For generating the minimax or Chebyshev polynomial coefficients required to approximate the custom trigonometric functions like $f_{BinBoot}$.

### 3.4 Custom Implementation Requirements (To Be Developed)

1. **Custom Parameter Set Generator**: A wrapper around `ckks.ParametersLiteral` that validates that the bottom-most moduli primes are properly ratio-locked to the scaling factor ($\Delta_0$), circumventing standard Lattigo defaults which usually enforce larger gaps.
2. **Polynomial Approximation Generators**: Go functions that output slice arrays of coefficients representing the Taylor/Chebyshev expansions of $f_{BinBoot}$ and $f_{G}$ (Table 4) bounded over the expected plaintext domains.
3. **BinBoot / GateBoot Evaluator Wrappers**:
* Implement a function `EvaluateBinBoot(eval ckks.Evaluator, ct *rlwe.Ciphertext)` that overrides the `EvalMod` step of Lattigo's default bootstrapper with the custom polynomial.
* Implement `EvaluateGateBoot(eval ckks.Evaluator, ct1, ct2 *rlwe.Ciphertext, gateType string)` which first computes `eval.Add(ct1, ct2, ctAdd)`, executes StC/CtS, extracts the real part, and evaluates the polynomial mapped to `gateType`.

4. **Accuracy Validation Module**: A test utility `ValidateBinaryTolerance(expected []float64, actual []float64)` that verifies all output floats are within the paper's specified error margin (e.g., testing that decrypted "1" logic gates are $> 0.75$ and "0" logic gates are $< 0.25$ before a final hard rounding step).