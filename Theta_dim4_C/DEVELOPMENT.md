# Theta_dim4 library in C

Work in progress. 

This library is a C implementation of the Sagemeath library Theta_dim4_sage to compute $2^e$-isogenies in dimension 4 (obtained by Kani's lemma) with level 2 theta coordinates based on the algorithms described in [1].

The C codebase is inspired from SQIsign NIST round 2 submission https://github.com/SQIsign/the-sqisign/tree/nist-v2. Unlike the latter, it is compiled with make and not makefile.

## Functionality

The goal of this library is to compute $2^e$-isogenies obtained with Kani's lemma of the form 

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
- sagemath to generate parameters for each choice of finite field.

As this library should be able to handle many finite fiels and to allow the user to generate new ones, modarith and sagemath were included for that purpose. These dependencies are unnecessary for applications when finite fields are already generated.

## How to run the code? [End goal]

### Step 1: parameter generation

The library should be able to generate and handle fields of the form $p=c\cdot 2^{e_2}\cdot \ell^{e_\ell}-1$, with $$c$$ small and $\ell$ prime. Note that $e_\ell$ can be zero (as in SQIsignHD).

Parameter generation scripts are located in the precomp_scripts file.

#### Generate a field to test Kani's 4-dimensional embedding

Go to `precomp_scripts` and type:

```sage generate_fp -p=<choice of prime> -name=<name of chosen prime> -test```

The prime `<choice of prime>` must be of the form $p=c\cdot 2^{e_2}\cdot \ell^{e_\ell}-1$ with $e_\ell>0$.

Then the script will generate:
- a file `src/gf/fp/fp_<name of chosen prime>.c` containing $\mathbb{F}_p$ arithmetic functions for the chosen prime `<choice of prime>`.
- a file `TODO` parameters $a_1, a_2, f_2, f_\ell$ such that $f_2\leq e_2-2$, $f_\ell\leq e_\ell$ and $a_1^2+a_2^2+\ell^{f_\ell}=2^{f_2}$.

#### Generate a field for other uses

Go to `precomp_scripts` and type:

```sage generate_fp -p=<prime> -name=<prime name>```

The chosen prime `<prime>` must be of the form $p=c\cdot 2^{e_2}\cdot \ell^{e_\ell}-1$ with $e_\ell=0$.

Then the script will generate a file `src/gf/fp/fp_<prime name>.c` containing $\mathbb{F}_p$ arithmetic functions for the chosen prime `<prime>`.

### Step 2: compilation and use

#### Testing Kani's 4-dimensional embedding

Go to `src` and type

```make test pname=<prime name> [options]```

To run the compiled code, type 

```./test_kani_4D_<prime name>```

This code should execute the equivalent `KaniEndo` and `KaniEndoHalf` functions (see `Theta_dim4_sage/README.md`) on the chosen base field indicated by `<prime name>`. The relevant parameters for `<prime name>` should already be generated (e.g. in Step 1).

TODO: define `[options]`.

#### Testing SIDH attacks

Go to `src` and type

```make wreck_SIKE pname=[p434, p503, p610, p751]```

To run the compiled code, type 

```./SIDH_attack_<prime name>```

## Structure of the library [with pending tasks]

The architecture of the library is greatly inspired from SQIsign-NIST-v2 https://github.com/SQIsign/the-sqisign/tree/nist-v2. The modules to be included in `src` are the following:

[Modules from SQIsign-NIST-v2 to include]

- `common`: code called in several files from SQIsign-NIST-v2, including pseudo-random number generators [directly taken from SQIsign-NIST-v2] [TODO].
- `gf`: finite field arithmetic ($\mathbb{F}_p$ and $\mathbb{F}_{p^2}$) [directly taken from SQIsign-NIST-v2] [IN PROGRESS].
- `ec`: elliptic curves along with their isogenies and pairings [directly taken from SQIsign-NIST-v2 with additional functions from src.old2 for odd degree isogenies and pairings] [TODO].
- `dim2`: for 2-dimensional isogenies [copy of `hd` in SQIsign-NIST-v2 with modifications in change of coordinate formulas] [TODO].
- `mp`: code for multiprecision arithmetic [directly taken from SQIsign-NIST-v2] [TODO].
- `params`: headers and C files with the parameters (depending on the chosen field) [Same role as `precomp` in SQIsign-NIST-v2] [TODO]. 

[Module(s) to create]

- `dim4`: for 4-dimensional theta structures and isogenies [implementation in C of functions from `Theta_dim4_sage/pkg/[isogenies,basis_change,theta_structures]`] [TODO].

The following main files are to be created in the `src` folder: 

- `kani_4D.c` (referenced in `kani_4D.h`) containing the equivalent of the functions `KaniEndo` and `KaniEndoHalf` from `Theta_dim4_sage/pkg/isogenies/Kani_endomorphism.py` [TODO].
- `test_kani_4D.c` (including `kani_4D.h`) testing the functions from `kani_4D.c` [TODO].
- `SIDH_attack.c` (including `kani_4D.h`) containing the SIDH attack code [TODO].






## References

[1] Pierrick Dartois, Fast computation of 2-isogenies in dimension 4 and cryptographic applications, In Journal of Algebra, Volume 683, December 2025. Preprint available at https://eprint.iacr.org/2024/1180.

[2] Pierrick Dartois, Antonin Leroux, Damien Robert and Benjamin Wesolowski, SQISignHD: New Dimensions in Cryptography, In Advances in Cryptology – EUROCRYPT 2024. Preprint available at https://eprint.iacr.org/2023/436.

[3] Pierrick Dartois, Jonathan Komada Eriksen, Tako Boris Fouotsa, Arthur Herlédan Le Merdy, Riccardo Invernizzi, Damien Robert, Ryan Rueger, Frederik Vercauteren and Benjamin Wesolowski, PEGASIS: Practical Effective Class Group Action using 4-Dimensional Isogenies, Accepted at CRYPTO 2025. Preprint available at https://eprint.iacr.org/2025/401.

[4] Unknown complete list of authors. qt-PEGASIS. Unpublished, 2025.