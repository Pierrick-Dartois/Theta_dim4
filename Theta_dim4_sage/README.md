# Theta_dim4 sagemath library

This library is a `python/sagemath` implementation of 4-dimensional 2-isogeny computations with the Theta model based on [2]. Only chains of 2-isogenies and Theta models of level 2 are implemented. 

Copyright (c) 2024, Pierrick Dartois.

The library is suited to the computation of isogenies between products of elliptic curves given by Kani's lemma, especially in the context of SIDH attacks and SQISignHD verification. We recover isogenies in dimension 1 given their evaluation on torsion points. Hence, the library can also be used to test SIDH attacks. Finally, this library also contains subroutines for effective group actions in the CSIDH context using 4-dimensional isogenies, as described in the PEGASIS paper [4]. Nonetheless, the PEGASIS API is neither in this repository nor documented here. We refer to https://github.com/pegasis4d/pegasis.


## Disclaimer

This library is experimental. Further tests and optimizations are still needed. 

Due to already implemented optimizations, the code is significantly different from what is described in Appendix F of the SQISignHD paper [1]. Another paper dedicated to dimension 4 $2$-isogeny computations in the Theta model describes more accurately the algorithms implemented in this library [2]. 

Note that this library heavily relies on dimension 2 isogeny computations in the Theta-model as described in [3].


## Requirement

You should have `python 3` and `sagemath` (version 10.0 at least) installed on your computer. Instructions to download `sagemath` can be found here https://doc.sagemath.org/html/en/installation/index.html.


## How to run the code?

Start a terminal and use the `cd` command to locate in the `Theta_dim4_sage` directory. From there, you can access two command line interfaces (CLI):
- `Tests.py` running generic tests on 4-dimensional 2-isogeny chains derived from Kani's lemma;
- `SIDH_attack.py` to test attacks against SIDH proposed in [5].

To use the CLI, type:

`sage <Tests.py/SIDH_attack.py> <arguments>`

More details on the `<arguments>` of each CLI are provided in the following. 


## <a name="Generic"></a> Generic tests

The script in `Tests.py` computes instanciations of 4-dimensional 2-isogeny chains derived from Kani's lemma given as endomorphisms:

$$F:=\left(\begin{matrix} a_1 & a_2 & \widehat{\sigma} & 0 \\ 
-a_2 & a_1 & 0 & \widehat{\sigma} \\ 
-\sigma & 0 & a_1 & -a_2 \\ 
0 & -\sigma & a_2 & a_1 \end{matrix}\right)\in \mathrm{End}(E_1^2\times E_2^2),$$

where:
- $E_1$ is a "random" supersingular elliptic curve, generated as the codomain of a random isogeny walk of degree $\simeq p$ starting from the supersingular elliptic curve $E_0$ of $j$-invariat $1728$;
- $\sigma: E_1\longrightarrow E_2$ is an isogeny of degree $\ell_B^{e_B}$ where $\ell_B$ is an odd prime ($\ell_B=3$ or $7$ here); 
- Integers $a_1, a_2, e_A$ have been chosen such that:
$$a_1^2+a_2^2+\ell_B^{e_B}=2^{e_A};$$
- The characteristic $p$ is of the form $p:=f\cdot 2^{f_A}\cdot \ell_B^{f_B}-1$, with $f_A\geq e_A+2$ and $f_B\geq e_B$ so that we have enough accessible torsion defined over $\mathbb{F}_{p^2}$ to compute $\sigma$ and evaluate it on $E_1[2^{e_A+2}]$.

This script includes two different ways of computing such a dimension 4 isogeny:
- the command line instruction `--KaniEndo` tests instances when we have the full torsion available.
- the function `--KaniEndoHalf` tests instances when only half the torsion is available (as in SQISignHD, see Section 4.3 of the paper).

### Parameters 

Different set of parameters can be found in the subdirectory `parameters`. The files labeled `parameters/parameters_3.txt` and `parameters/parameters_7.txt` contain lists of parameters for the values $\ell_B=3$ and $\ell_B=7$ respectively. Every line of those files contains a list:

$$[e_A, e_B, a_1, a_2, f, f_A, f_B, p(, m)]$$

such that:
- $p:=f\cdot 2^{f_A}\cdot \ell_B^{f_B}-1$;
- $f_A\geq e_A+2$;
- $f_B\geq e_B$;
- $m=\max(v_2(a_1),v_2(a_2))$ (this parameter is not present in the list when $\ell_B=3$ since it is always $1$ in this case).

These parameters have been generated with the functions `find_param` and `find_param_gen` of `parameters/parameter_generation.py`. 

#### Displaying parameters 

You can display available parameters with the command:

```
sage Tests.py --display
```

or in short:

```
sage Tests.py -d
```

This command will display all saved parameters with their index in the list. Those indices will be used as reference to execute tests on a specific set of parameters.

If you want to restrict the displayed lists of parameters, type:

```
sage Tests.py -d -l_B=<value>
```

to display all parameters with a given `<value>` of $\ell_B$ or type:

```
sage Tests.py -d -l_B=<value l_B> -i=<value index>
```

to display the list of parameters with a given `<value l_B>` of $\ell_B$ indexed with `<value index>`.

#### Creating new parameters

If you want to create new lists of parameters, you can type:

```sage --add_params -l_B=<value l_B> -e_A=<value e_A>```

this command will create lists of parameters with $\ell_B=$`<value l_B>`, $e_A=$`<value e_A>` and $e_B$ varying in the range $2^{e_A/2}\leq \ell_B^{e_B}\leq 2^{e_A}$. Depending on the specified values, multiple parameters or no parameter can be found. The new parameters are saved in an existing file (after the previously saved parameters) or in a newly created file `parameters/parameters_<value l_B>.txt`. The newly created parameters are displayed along with their index in the file.

EXAMPLE:
```
% sage Tests.py --add_params -l_B=11 -e_A=66
===========================================
New parameters with second prime l_B=11.
===========================================

 - Index in the list of parameters = 0
 - Prime characteristic p = 1 * 2**70 * 11**11 - 1
 - Degree of the embedded isogeny sigma q = 11**11
 - a1 = 1440875693
 - a2 = -8468226098
 - m = max(v_2(a1),v_2(a2)) = 1
 - Length of the 4-dimensional 2-isogeny = 66
 ```

The parameter generation is subexponential and may take a while. We do not recommend to select `<value e_A>` bigger than $200$.

### Chain with the full torsion available

When we can access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the whole $2$-isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can be computed directly. The command line instruction `--KaniEndo` tests this direct chain computation for specified sets of parameters. 

Type:

```
sage Tests.py --KaniEndo
```

to test all saved lists of parameters. 

Type:

```
sage Tests.py --KaniEndo -l_B=<value l_B>
```

to test all saved lists of parameters with $\ell_B=$`<value l_B>`.

Type:

```
sage Tests.py --KaniEndo -l_B=<value l_B> -i=<value index>
```

to test the list of parameters with $\ell_B=$`<value l_B>` indexed by `<value index>`.

#### For experimented users only

You can test a specific set of parameters chosen manually. If the set of parameters is not well chosen, this could lead to errors.

```
sage Tests.py --KaniEndo -l_B=<value l_B> -e_A=<value e_A> -e_B=<value e_B> -a1=<value a1> -a2=<value a2> -f=<value f> -f_A=<value f_A> -f_B=<value f_B> -p=<value p> [optional: -m=<value m>]
```

EXAMPLE:
```
% sage Tests.py --KaniEndo -l_B=7 -e_A=17 -e_B=3 -a1=123 -a2=-340 -f=3 -f_A=20 -f_B=3 -p=1078984703
Testing KaniEndo with parameters:
 - Prime characteristic p = 3 * 2**20 * 7**3 - 1
 - Degree of the embedded isogeny sigma q = 7**3
 - a1 = 123
 - a2 = -340
 - m = max(v_2(a1),v_2(a2)) = 2
 - Length of the dimension 4 2-isogeny = 17
Setup: 0.04501986503601074 s
Random walk: 0.015614986419677734 s
Generation of sigma: 0.004041910171508789 s
Generation and evaluation of the torsion basis: 0.002375364303588867 s
Strategy computation: 0.0004138946533203125 s
Dimension 4 endomorphism: 0.16916584968566895 s
Is evaluation correct?
True
```

### Chain with half the necessary torsion (as in SQISignHD)

When we cannot access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the computation of the 2-isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can has to be divided in two, as in Section 4.3 of the SQISignHD paper. Namely, we compute two isogeny chains $F_1: E_1^2\times E_2^2\longrightarrow C$ and $\widetilde{F_2}: E1^2\times E2^2\longrightarrow C$ such that $F=F_2\circ F_1$. The command line instruction `--KaniEndoHalf` tests this computation for specified sets of parameters. 

Type:

```
sage Tests.py --KaniEndoHalf
```

to test all saved lists of parameters. 

Type:

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B>
```

to test all saved lists of parameters with $\ell_B=$`<value l_B>`.

Type:

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B> -i=<value index>
```

to test the list of parameters with $\ell_B=$`<value l_B>` indexed by `<value index>`.

#### For experimented users only

You can test a specific set of parameters chosen manually. If the set of parameters is not well chosen, this could lead to errors.

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B> -e_A=<value e_A> -e_B=<value e_B> -a1=<value a1> -a2=<value a2> -f=<value f> -f_A=<value f_A> -f_B=<value f_B> -p=<value p> [optional: -m=<value m>]
```

EXAMPLE:
```
% sage Tests.py --KaniEndoHalf -l_B=7 -e_A=17 -e_B=3 -a1=123 -a2=-340 -f=3 -f_A=20 -f_B=3 -p=1078984703
Testing KaniEndoHalf with parameters:
 - Prime characteristic p = 3 * 2**20 * 7**3 - 1
 - Degree of the embedded isogeny sigma q = 7**3
 - a1 = 123
 - a2 = -340
 - m = max(v_2(a1),v_2(a2)) = 2
 - Length of the dimension 4 2-isogeny = 17
 - Used available torsion = 2**11
Setup: 0.04306387901306152 s
Random walk: 0.015376091003417969 s
Generation of sigma: 0.004205942153930664 s
Generation and evaluation of the torsion basis: 0.002997159957885742 s
Computation of strategies: 0.00025773048400878906 s
Dimension 4 endomorphism: 0.20006227493286133 s
Is evaluation correct?
True
Time evaluation: 0.00693202018737793 s
```

### Implementation restrictions

When computing dimension 4 isogenies with the Theta model, extra care is needed to compute gluing isogenies (isogenies starting from a product of abelian varieties of dimension less than 4). In most cases, when we compute an isogeny chain derived from Kani's lemma, we only encounter gluings at the beginning. However, it may happen in small characteristic that an isogeny splits into a product of abelian varieties in the middle of the chain. In this case, we have to compute a gluing afterwards and this step generally fails since we can only compute gluings when we know their location in advance in the chain. Our code then returns a NotImplementedError.    

In particular, the set of parameters with $\ell_B=3$ and index 0 always strikes a failure.

```
% sage Tests.py --KaniEndo -l_B=3 -i=0
Testing KaniEndo with parameters:
 - Prime characteristic p = 1 * 2**18 * 3**5 - 1
 - Degree of the embedded isogeny sigma q = 3**5
 - a1 = 238
 - a2 = -93
 - m = max(v_2(a1),v_2(a2)) = 1
 - Length of the dimension 4 2-isogeny = 16
Setup: 0.04303693771362305 s
Random walk: 0.018531084060668945 s
Generation of sigma: 0.004489898681640625 s
Generation and evaluation of the torsion basis: 0.002763986587524414 s
Strategy computation: 0.0003800392150878906 s
[...]
NotImplementedError: The codomain of this 2-isogeny could not be computed.
We may have encountered a product of abelian varieties
somewhere unexpected along the chain.
This is exceptionnal and should not happen in larger characteristic.
```


## SIDH attacks

The file `SIDH_attack.py` runs a full key recovery attack against the Supersingular Isogeny Diffie Hellman (SIDH) key exchange following the approach of [2, § 5.5] for all SIKE NIST primes (p434, p503, p610, p751) and with an arbitrary starting curve $E_1$.

### What does this script do?

For a selected SIKE NIST prime (p434, p503, p610, p751) of the form $p=2^{e_2}\cdot 2^{e_3}-1$, the scripts samples a supersingular elliptic curve $E_1$ defined over $\mathbb{F}_{p^2}$ (starting from the curve of $j$-invariant $1728$ and walking a random isogeny path of degree $\approx p$). Then it performs an SIDH key exchange starting from $E_1$ and runs an attack to recover the shared secret key $E_{AB}$ of Alice and Bob.

#### The SIDH protocol:

- Let $(P_A,Q_A)$ and $(P_B,Q_B)$ be basis of $E_1[2^{e_2}]$ and $E_1[3^{e_3}]$ respectively. 
- Alice samples $s_A\in\{0,\cdots, 2^{e_2}-1\}$ and computes $\varphi_A: E_1\longrightarrow E_A$ of kernel $\langle P_A+[s_A]Q_A\rangle$ and sends $(\varphi_A(P_B),\varphi_A(Q_B))$ to Bob. 
- Bob samples $s_B\in\{0,\cdots, 3^{e_3}-1\}$ and computes $\varphi_B: E_1\longrightarrow E_B$ of kernel $\langle P_B+[s_B]Q_B\rangle$ and sends $(\varphi_B(P_A),\varphi_B(Q_A))$ to Alice.
- Alice computes $\psi_A: E_B\longrightarrow E_{BA}$ of kernel $\langle \varphi_B(P_A)+[s_A]\varphi_B(Q_A)\rangle$.
- Bob computes $\psi_B: E_A\longrightarrow E_{AB}\simeq E_{BA}$ of kernel $\langle \varphi_A(P_B)+[s_B]\varphi_A(Q_B)\rangle$.


#### The key recovery attack:

The attacker uses $(\varphi_B(P_A),\varphi_B(Q_A))$ to compute the secret $s_B$.
- They use precomputed integers $e, a_1, a_2$ such that $a_1^2+a_2^2+3^{e_3}=2^e$ (stored in `SIDH/parameters_SIKE_NIST.txt`). $e>e_2-2$ so the available torsion is insufficient to use `KaniEndo` but $\lceil e/2\rceil+2\leq e_2$ so `KaniEndoHalf` can (and has to) be used.
- Using `KaniEndoHalf` and $(\varphi_B(P_A),\varphi_B(Q_A))$, the attacker can evaluate $\varphi_B$ everywhere and especially on $(P_B,Q_B)$.
- With discrete logarithms, they can recover $\ker(\varphi_B)$ so $s_B$.
- They can then compute $\psi_B: E_A\longrightarrow E_{AB}$ of kernel $\langle \varphi_A(P_B)+[s_B]\varphi_A(Q_B)\rangle$ and obtain the shared secret key $E_{AB}$.

### Instructions to run the script

To run the complete attack, type:

```
sage SIDH_attack.py --Protocol_and_Attack -p=prime
```

where `prime` is a choice of SIKE NIST prime (`p434`, `p503`, `p610` or `p751`).

To display the available primes, type:

```
sage SIDH_attack.py --display
```

or in short format:

```
sage SIDH_attack.py -d
```

If you want to run an SIDH key exchange without the attack, type:

```
sage SIDH_attack.py --Protocol -p=prime
```

where `prime` is again a choice of SIKE NIST prime (`p434`, `p503`, `p610` or `p751`).


## Benchmarks

The user can reproduce the timings from [2, Tables 2-3] comparing 4-dimensional 2-isogeny chains computations (with `KaniEndo` and `KaniEndoHalf`) and 1-dimensional 2-isogeny chains of the same length for some parameters contained in the `parameters` file. Simply type:

```
sage Benchmarks.py
```

The results are then saved in the file `Benchmarking_results.csv`. Note that for each parameter, tests are reproduced 100 times so the benchmarking may take more than 1h.

The user can also benchmark the SIDH attacks and reproduce [2, Table 4]. Simply type:

```
sage Benchmarks_SIDH.py
```

The results are then saved in the file `Benchmarking_results_SIDH.csv`. Note that for each parameter, tests are reproduced 100 times so the benchmarking may take several tens of minutes.


## Organization of the library

The main test files `Tests.py`, `Benchmarks.py`, `Tests.py`, `SIDH_attack.py` and `Benchmarks_SIDH.py` are located in the main directory `Theta_dim4_sage` of this library. 

The main test files all imports functions from a package `pkg` containing several modules:
- `basis_change` contains code for changing level 2 Theta structures in dimension 2 and 4 useful for the computation of gluing and splitting isogenies (in the beginning and in the end of the chain).
- `isogenies` contains code to compute 2-isogenies and chains of 2-isogenies in dimension 4 in the Theta model.
- `isogenies_dim2` contains code to compute 2-isogenies and chains of 2-isogenies in dimension 2 in the Theta model, based on the implementation of [3]. This code is useful for the first steps of dimension 4 2-isogeny chains (involving gluings).
- `montgomery_isogenies` contains code to compute dimension 1 isogeny on the Kummer line faster than `sagemath` using $x$-only arithmetic. These files are due to Giacomo Pope.
- `parameters` contains the parameters mentionned in [Generic tests](#Generic).
- `theta_structures` contains code for different models in dimension 1, 2 and 4 and translations between those models (Montgomery model on elliptic curves (dimension 1) or product of elliptic curves, level 2 Theta models in dimensions 1, 2 and 4).
- `SIDH` contains data and parameters to run the SIDH attacks.
- `utilities` contains several useful functions to work on supersingular elliptic curves faster than with the `sagemath` generic functions (discrete logarithms, Weil pairings...).

## License

`Theta_dim4` is licensed under Apache-2.0. See LICENSE and NOTICE in the root directory.

Third party code is used in this directory (`Theta_dim4_sage`):

- `montgomery_isogenies`; MIT: "Copyright (c) 2023 Giacomo Pope"
- `utilities`; "Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope"
- `isogenies_dim2`; MIT: "Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert" 

## References

[1] Pierrick Dartois, Antonin Leroux, Damien Robert and Benjamin Wesolowski, SQISignHD: New Dimensions in Cryptography, In Advances in Cryptology – EUROCRYPT 2024. https://eprint.iacr.org/2023/436.

[2] Pierrick Dartois, Fast computation of 2-isogenies in dimension 4 and cryptographic applications, IACR Cryptology ePrint Archive, Paper 2024/1180, 2024. https://eprint.iacr.org/2024/1180

[3] Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert, An Algorithmic Approach to (2,2)-isogenies in the Theta Model and Applications to Isogeny-based Cryptography, In Advances in Cryptology – ASIACRYPT 2024. https://eprint.iacr.org/2023/1747.

[4] Pierrick Dartois, Jonathan Komada Eriksen, Tako Boris Fouotsa, Arthur Herlédan Le Merdy, Riccardo Invernizzi, Damien Robert, Ryan Rueger, Frederik Vercauteren and Benjamin Wesolowski, PEGASIS: Practical Effective Class Group Action using 4-Dimensional Isogenies, IACR Cryptology ePrint Archive, Paper 2025/401, 2025. https://eprint.iacr.org/2025/401.

[5] Damien Robert. Breaking SIDH in polynomial time, In Advances in Cryptology – EUROCRYPT 2023. https://eprint.iacr.org/2022/1038

