# Theta_dim4 library in C

Work in progress. 

This library is a C implementation of the Sagemeath library Theta_dim4_sage to compute $2^e$-isogenies in dimension 4 (obtained by Kani's lemma) with level 2 theta coordinates based on the algorithms described in [1].

The C codebase is inspired from SQIsign NIST round 2 submission https://github.com/SQIsign/the-sqisign/tree/nist-v2. Unlike the latter, it is compiled with make and not makefile.

## Functionality

The goal of this library is to compute $2^e$-isogenies of the form 

$$F:=\left(\begin{matrix} a_1 & a_2 & \widehat{\sigma} & 0 \\ 
-a_2 & a_1 & 0 & \widehat{\sigma} \\ 
-\sigma & 0 & a_1 & -a_2 \\ 
0 & -\sigma & a_2 & a_1 \end{matrix}\right)\in \mathrm{End}(E_1^2\times E_2^2),$$

where $\sigma: E_1\longrightarrow E_2$ is an isogeny and $a_1, a_2\in\mathbb{Z}$ are such that $a_1^2+a_2^2+\deg(\sigma)=2^e$. 

Such isogenies can be used to interpolate $\sigma$ given its image on $2^*$-torsion points. This library:
- tests this functionality for different base fields to reproduce the timings in [1, Tables 2-3].
- applies this functionality to SIDH attacks at all security levels.
This can also be directly used directly in SQIsignHD [2] verification (with a submodule inclusion).

A minimal copy of this library (with less dependencies and finite fields) could also be easily adapted to implement (qt-)PEGASIS [3,4].

## Dependencies

- C11-compatible compiler.
- modarith (https://github.com/mcarrickscott/modarith), a library to automatically generate code for the $\mathbb{F}_p$ arithmetic.

As this library should be able to handle many finite fiels and to allow the user to generate new ones, modarith was included for that purpose.

## How to run the code? End goal

### Step 1: parameter generation



### Step 2: compilation and use

#### Testing Kani's 4D embedding

#### Testing SIDH attacks





## References

[1] Pierrick Dartois, Fast computation of 2-isogenies in dimension 4 and cryptographic applications, In Journal of Algebra, Volume 683, December 2025. Preprint available at https://eprint.iacr.org/2024/1180.

[2] Pierrick Dartois, Antonin Leroux, Damien Robert and Benjamin Wesolowski, SQISignHD: New Dimensions in Cryptography, In Advances in Cryptology – EUROCRYPT 2024. Preprint available at https://eprint.iacr.org/2023/436.

[3] Pierrick Dartois, Jonathan Komada Eriksen, Tako Boris Fouotsa, Arthur Herlédan Le Merdy, Riccardo Invernizzi, Damien Robert, Ryan Rueger, Frederik Vercauteren and Benjamin Wesolowski, PEGASIS: Practical Effective Class Group Action using 4-Dimensional Isogenies, Accepted at CRYPTO 2025. Preprint available at https://eprint.iacr.org/2025/401.

[4] Unknown complete list of authors. qt-PEGASIS. Unpublished, 2025.