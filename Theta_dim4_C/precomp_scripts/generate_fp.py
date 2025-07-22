#!/usr/bin/env python3

from sage.functions.other import ceil
from sage.arith.misc import factor
from sage.arith.misc import inverse_mod
from sage.rings.factorint import factor_trial_division
from sage.structure.factorization import Factorization
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import Mod

import os
import argparse
from shutil import which
from textwrap import dedent
from re import match


def prime_shape(p):
    """Factor p + 1 = c * 2^e_2 * ell^e_ell, returning tuple (c, e_2, ell, e_ell)

    - `c` is a (usually small) cofactor
    - `ell` the odd prime dividing p + 1 with the highest multiplicity (`e_ell`)
    """
    assert (p + 1) % 2 == 0

    # We choose `ell` by finding the prime with the highest multiplicity
    # First sort by prime
    factors = sorted(factor(p + 1), key=lambda x: x[0])

    # Get multiplicity of 2 in power (`factor` returns a sorted list)
    e_2 = factors[0][1]

    # Sort remaining by multiplicity (ignoring 2)
    factors = sorted(factors[1:], key=lambda x: x[1], reverse=True)

    if factors:
        ell, e_ell = factors[0]
        c = (p + 1) / (2**e_2 * ell**e_ell)
    else:
        # p is a mersenne prime
        ell, e_ell = 1, 1
        c = 1

    assert p + 1 == c * 2**e_2 * ell**e_ell

    return c, e_2, ell, e_ell


def is_sum_of_two_squares(N, smoothness=1000):
    r"""Return factorisation of N if it is a sum of two squares and `smoothness`-smooth, False otherwise"""

    partial_factorization = factor_trial_division(N, smoothness)

    for p, e_p in partial_factorization:
        if not (p == 2 or p % 4 == 1 or (p % 4 == 3 and e_p % 2 == 0)):
            return False

    return partial_factorization


# Taken from PEGASIS source
# https://github.com/pegasis4d/pegasis
def two_squares_factored(factors):
    """
    This is the function `two_squares` from sage, except we give it the
    factorisation of n already.
    """
    F = Factorization(factors)
    for p, e in F:
        if e % 2 == 1 and p % 4 == 3:
            raise ValueError("%s is not a sum of 2 squares" % n)

    n = F.expand()
    if n == 0:
        return (0, 0)
    a = ZZ.one()
    b = ZZ.zero()
    for p, e in F:
        if p == 1:
            continue
        if e >= 2:
            m = p ** (e // 2)
            a *= m
            b *= m
        if e % 2 == 1:
            if p == 2:
                # (a + bi) *= (1 + I)
                a, b = a - b, a + b
            else:  # p = 1 mod 4
                # Find a square root of -1 mod p.
                # If y is a non-square, then y^((p-1)/4) is a square root of -1.
                y = Mod(2, p)
                while True:
                    s = y ** ((p - 1) / 4)
                    if not s * s + 1:
                        s = s.lift()
                        break
                    y += 1
                # Apply Cornacchia's algorithm to write p as r^2 + s^2.
                r = p
                while s * s > p:
                    r, s = s, r % s
                r %= s

            # Multiply (a + bI) by (r + sI)
            a, b = a * r - b * s, b * r + a * s

    a = a.abs()
    b = b.abs()
    assert a * a + b * b == n
    if a <= b:
        return (a, b)
    else:
        return (b, a)


def find_embedding_params(p_shape):
    r"""
    Given exponents e2_max and el_max of the accessible 2**e2_max and
    l**el_max-torsion respectively, this function finds integers e2, el, a1, a2
    such that 2**e2=l**el+a1**2+a2**2, where el is the biggest possible.
    """
    _, e2_max, l, el_max = p_shape

    e2 = e2_max
    el = el_max

    for el in range(el_max, 0, -1):
        for e2 in range(e2_max):
            N = 2**e2 - l**el
            if N > 0:
                if N_factored := is_sum_of_two_squares(N):
                    a1, a2 = two_squares_factored(N_factored)
                    return e2, el, a1, a2

    return None


def int_to_montgemery_fp_const(x, p, Nlimbs, Radix):
    R = 2 ** (Nlimbs * Radix)
    el = (x * R) % p
    vs = [(int(el) >> Radix * i) % 2**Radix for i in range(Nlimbs)]
    return "{" + ", ".join(map(hex, vs)) + "}"


def write_field_file(p):
    d_word_params = {}

    params = ["Wordlength", "Nlimbs", "Radix", "Nbits", "Nbytes"]

    for line in open("field.c").readlines():
        for param in params:
            if m := match(f"^#define {param} ([0-9]*)$", line):
                d_word_params[param] = int(m.group(1))

    Nlimbs = d_word_params["Nlimbs"]
    Radix = d_word_params["Radix"]

    lines = []

    lines += [""]
    lines += ["/******************************************************************************"]
    lines += ["API functions calling generated code above"]
    lines += ["******************************************************************************/"]
    lines += [""]

    lines += [f"const digit_t ZERO[NWORDS_FIELD] = {int_to_montgemery_fp_const(0, p, Nlimbs, Radix)};"]
    lines += [f"const digit_t ONE[NWORDS_FIELD] = {int_to_montgemery_fp_const(1, p, Nlimbs, Radix)};"]
    lines += [
        f"static const digit_t TWO_INV[NWORDS_FIELD] = {int_to_montgemery_fp_const(inverse_mod(2, p), p, Nlimbs, Radix)};"
    ]
    lines += [
        f"\nstatic const digit_t THREE_INV[NWORDS_FIELD] = {int_to_montgemery_fp_const(inverse_mod(3, p), p, Nlimbs, Radix)};"
    ]

    lines += [
        dedent(
            """
        void
        fp_set_small(fp_t *x, const digit_t val)
        {
            modint((int)val, *x);
        }

        void
        fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)
        {
            modmli(*a, (int)val, *x);
        }

        void
        fp_set_zero(fp_t *x)
        {
            modzer(*x);
        }

        void
        fp_set_one(fp_t *x)
        {
            modone(*x);
        }

        uint32_t
        fp_is_equal(const fp_t *a, const fp_t *b)
        {
            return -(uint32_t)modcmp(*a, *b);
        }

        uint32_t
        fp_is_zero(const fp_t *a)
        {
            return -(uint32_t)modis0(*a);
        }

        void
        fp_copy(fp_t *out, const fp_t *a)
        {
            modcpy(*a, *out);
        }

        void
        fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
        {
            modcsw((int)(ctl & 0x1), *a, *b);
        }

        void
        fp_add(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modadd(*a, *b, *out);
        }

        void
        fp_sub(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modsub(*a, *b, *out);
        }

        void
        fp_neg(fp_t *out, const fp_t *a)
        {
            modneg(*a, *out);
        }

        void
        fp_sqr(fp_t *out, const fp_t *a)
        {
            modsqr(*a, *out);
        }

        void
        fp_mul(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modmul(*a, *b, *out);
        }

        void
        fp_inv(fp_t *x)
        {
            modinv(*x, NULL, *x);
        }

        uint32_t
        fp_is_square(const fp_t *a)
        {
            return -(uint32_t)modqr(NULL, *a);
        }

        void
        fp_sqrt(fp_t *a)
        {
            modsqrt(*a, NULL, *a);
        }

        void
        fp_half(fp_t *out, const fp_t *a)
        {
            modmul(TWO_INV, *a, *out);
        }

        void
        fp_exp3div4(fp_t *out, const fp_t *a)
        {
            modpro(*a, *out);
        }

        void
        fp_div3(fp_t *out, const fp_t *a)
        {
            modmul(THREE_INV, *a, *out);
        }

        void
        fp_encode(void *dst, const fp_t *a)
        {
            modexp(a, dst);
        }

        uint32_t
        fp_decode(fp_t *d, const void *src)
        {
            return modimp(src,d);
        }
    """
        )
    ]

    with open("field.c", "a+") as file:
        file.writelines([line + "\n" for line in lines])

    return d_word_params


def factor_pp1(p, factor_l=True):
    L = list(factor(p + 1))
    e2 = L[0][1]

    if factor_l:
        il = 1
        for i in range(1, len(L)):
            if L[i][1] > L[il][1]:
                il = i
        l, el = L[il]

        c = (p + 1) // (2**e2 * l**el)

        return c, e2, l, el
    else:
        c = (p + 1) // (2**e2)

        return c, e2


class Ibz:
    def __init__(self, v):
        self.v = int(v)

    def _literal(self, sz):
        val = int(self.v)
        sgn = val < 0
        num_limbs = (abs(val).bit_length() + sz - 1) // sz if val else 0
        limbs = [(abs(val) >> sz * i) & (2**sz - 1) for i in range(num_limbs or 1)]
        data = {
            "._mp_alloc": 0,
            "._mp_size": (-1) ** sgn * num_limbs,
            "._mp_d": "(mp_limb_t[]) {" + ",".join(map(hex, limbs)) + "}",
        }
        return "{{" + ", ".join(f"{k} = {v}" for k, v in data.items()) + "}}"


def write_constants_file(p, name, d_word_params, p_shape, args=None):
    filename = f"../src/params/fp_params/constants_{name}.c"

    c, e_2, ell, e_ell = p_shape

    lines = []
    lines += ["#include <constants.h>"]
    lines += [""]
    lines += [f"const uint64_t NWORDS_FIELD = {d_word_params["Radix"]};"]
    lines += [f"const uint64_t NWORDS_ORDER = {ceil(d_word_params["Nbits"] / d_word_params["Wordlength"])};"]

    charac = Ibz(p)
    tor_even = Ibz(2**e_2)
    tor_odd = Ibz((p + 1) // 2**e_2)

    lines += [f"const uint64_t TORSION_EVEN_POWER = {e_2};"]
    lines += [f"const ibz_t CHARACTERISTIC = {charac._literal(64)};"]
    lines += [f"const ibz_t TORSION_EVEN = {tor_even._literal(64)};"]
    lines += [f"const ibz_t TORSION_ODD = {tor_odd._literal(64)};"]

    if args is not None:
        f_2, f_l, a1, a2 = args

        tor_l = Ibz(ell**e_ell)
        comp_l = Ibz(c * 2**e_2)
        A1 = Ibz(a1)
        A2 = Ibz(a2)

        lines += [f"const uint64_t L = {ell};"]
        lines += [f"const uint64_t TORSION_L_POWER = {e_ell};"]
        lines += [f"const ibz_t TORSION_L = {tor_l._literal(64)};"]
        lines += [f"const ibz_t TORSION_COMPLEMENT_L = {comp_l._literal(64)};"]
        lines += [f"const uint64_t CONST_FL = {f_l};"]
        lines += [f"const uint64_t CONST_F2 = {f_2};"]
        lines += [f"const ibz_t CONST_A1 = {A1._literal(64)};"]
        lines += [f"const ibz_t CONST_A2 = {A2._literal(64)};"]

    print(f"Writing to file: {filename}")
    with open(filename, "w") as file:
        file.writelines([line + "\n" for line in lines])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--prime")
    parser.add_argument("-n", "--name")
    parser.add_argument("-t", "--test", action="store_true")

    args = parser.parse_args()

    # `addchain` must be compiled, and so cannot be shipped as-is
    if not which("addchain"):
        print("Cannot find `addchain` in `$PATH`. This is necessary for `modarith`")
        exit(1)

    # Bad habit to use `eval`, but since data is not supplied by untrused user it is probably okay
    # Could be replaced with https://stackoverflow.com/a/69540962
    p = eval(args.prime)
    p_shape = prime_shape(p)

    addl_args = None

    if args.test:
        addl_args = find_embedding_params(p_shape)

    os.system(f"python ../external/modarith/monty.py 64 {p} > /dev/null 2>&1")

    d_word_params = write_field_file(p)
    write_constants_file(p, args.name, d_word_params, p_shape, args=addl_args)

    os.system("mv -v field.c ../src/gf/fp/fp_" + args.name + ".c")
    os.system("rm -v time.c")
