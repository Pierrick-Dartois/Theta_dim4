
// Automatically generated modular arithmetic C code
// Command line : python ../modarith/monty.py 64
// 6763808680755523643471605336644762692711582514005431939938907264901915992959672694270019729488822011354742259711
// Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)

#include <stdint.h>
#include <stdio.h>

#define sspint int64_t
#define spint uint64_t
#define udpint __uint128_t
#define dpint __uint128_t

#define Wordlength 64
#define Nlimbs 7
#define Radix 54
#define Nbits 372
#define Nbytes 47

#define MONTGOMERY
// propagate carries
static spint inline prop(spint *n) {
  int i;
  spint mask = ((spint)1 << 54u) - (spint)1;
  sspint carry = (sspint)n[0];
  carry >>= 54u;
  n[0] &= mask;
  for (i = 1; i < 6; i++) {
    carry += (sspint)n[i];
    n[i] = (spint)carry & mask;
    carry >>= 54u;
  }
  n[6] += (spint)carry;
  return -((n[6] >> 1) >> 62u);
}

// propagate carries and add p if negative, propagate carries again
static int inline flatten(spint *n) {
  spint carry = prop(n);
  n[0] -= (spint)1u & carry;
  n[3] += ((spint)0x2488d500000000u) & carry;
  n[4] += ((spint)0x25cfee16507a30u) & carry;
  n[5] += ((spint)0x22c7bc80527adeu) & carry;
  n[6] += ((spint)0xb4000a89caeau) & carry;
  (void)prop(n);
  return (int)(carry & 1);
}

// Montgomery final subtract
static int inline modfsb(spint *n) {
  n[0] += (spint)1u;
  n[3] -= (spint)0x2488d500000000u;
  n[4] -= (spint)0x25cfee16507a30u;
  n[5] -= (spint)0x22c7bc80527adeu;
  n[6] -= (spint)0xb4000a89caeau;
  return flatten(n);
}

// Modular addition - reduce less than 2p
static void inline modadd(const spint *a, const spint *b, spint *n) {
  spint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[5] = a[5] + b[5];
  n[6] = a[6] + b[6];
  n[0] += (spint)2u;
  n[3] -= (spint)0x4911aa00000000u;
  n[4] -= (spint)0x4b9fdc2ca0f460u;
  n[5] -= (spint)0x458f7900a4f5bcu;
  n[6] -= (spint)0x16800151395d4u;
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[3] += ((spint)0x4911aa00000000u) & carry;
  n[4] += ((spint)0x4b9fdc2ca0f460u) & carry;
  n[5] += ((spint)0x458f7900a4f5bcu) & carry;
  n[6] += ((spint)0x16800151395d4u) & carry;
  (void)prop(n);
}

// Modular subtraction - reduce less than 2p
static void inline modsub(const spint *a, const spint *b, spint *n) {
  spint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  n[5] = a[5] - b[5];
  n[6] = a[6] - b[6];
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[3] += ((spint)0x4911aa00000000u) & carry;
  n[4] += ((spint)0x4b9fdc2ca0f460u) & carry;
  n[5] += ((spint)0x458f7900a4f5bcu) & carry;
  n[6] += ((spint)0x16800151395d4u) & carry;
  (void)prop(n);
}

// Modular negation
static void inline modneg(const spint *b, spint *n) {
  spint carry;
  n[0] = (spint)0 - b[0];
  n[1] = (spint)0 - b[1];
  n[2] = (spint)0 - b[2];
  n[3] = (spint)0 - b[3];
  n[4] = (spint)0 - b[4];
  n[5] = (spint)0 - b[5];
  n[6] = (spint)0 - b[6];
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[3] += ((spint)0x4911aa00000000u) & carry;
  n[4] += ((spint)0x4b9fdc2ca0f460u) & carry;
  n[5] += ((spint)0x458f7900a4f5bcu) & carry;
  n[6] += ((spint)0x16800151395d4u) & carry;
  (void)prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 2828534504801797697349267301285903
// Modular multiplication, c=a*b mod 2p
static void inline modmul(const spint *a, const spint *b, spint *c) {
  dpint t = 0;
  spint p3 = 0x2488d500000000u;
  spint p4 = 0x25cfee16507a30u;
  spint p5 = 0x22c7bc80527adeu;
  spint p6 = 0xb4000a89caeau;
  spint q = ((spint)1 << 54u); // q is unsaturated radix
  spint mask = (spint)(q - (spint)1);
  t += (dpint)a[0] * b[0];
  spint v0 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  spint v1 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  spint v2 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  t += (dpint)v0 * (dpint)p3;
  spint v3 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  t += (dpint)v0 * (dpint)p4;
  t += (dpint)v1 * (dpint)p3;
  spint v4 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[5];
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)a[5] * b[0];
  t += (dpint)v0 * (dpint)p5;
  t += (dpint)v1 * (dpint)p4;
  t += (dpint)v2 * (dpint)p3;
  spint v5 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[0] * b[6];
  t += (dpint)a[1] * b[5];
  t += (dpint)a[2] * b[4];
  t += (dpint)a[3] * b[3];
  t += (dpint)a[4] * b[2];
  t += (dpint)a[5] * b[1];
  t += (dpint)a[6] * b[0];
  t += (dpint)v0 * (dpint)p6;
  t += (dpint)v1 * (dpint)p5;
  t += (dpint)v2 * (dpint)p4;
  t += (dpint)v3 * (dpint)p3;
  spint v6 = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[1] * b[6];
  t += (dpint)a[2] * b[5];
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)a[5] * b[2];
  t += (dpint)a[6] * b[1];
  t += (dpint)v1 * (dpint)p6;
  t += (dpint)v2 * (dpint)p5;
  t += (dpint)v3 * (dpint)p4;
  t += (dpint)v4 * (dpint)p3;
  c[0] = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[2] * b[6];
  t += (dpint)a[3] * b[5];
  t += (dpint)a[4] * b[4];
  t += (dpint)a[5] * b[3];
  t += (dpint)a[6] * b[2];
  t += (dpint)v2 * (dpint)p6;
  t += (dpint)v3 * (dpint)p5;
  t += (dpint)v4 * (dpint)p4;
  t += (dpint)v5 * (dpint)p3;
  c[1] = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[3] * b[6];
  t += (dpint)a[4] * b[5];
  t += (dpint)a[5] * b[4];
  t += (dpint)a[6] * b[3];
  t += (dpint)v3 * (dpint)p6;
  t += (dpint)v4 * (dpint)p5;
  t += (dpint)v5 * (dpint)p4;
  t += (dpint)v6 * (dpint)p3;
  c[2] = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[4] * b[6];
  t += (dpint)a[5] * b[5];
  t += (dpint)a[6] * b[4];
  t += (dpint)v4 * (dpint)p6;
  t += (dpint)v5 * (dpint)p5;
  t += (dpint)v6 * (dpint)p4;
  c[3] = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[5] * b[6];
  t += (dpint)a[6] * b[5];
  t += (dpint)v5 * (dpint)p6;
  t += (dpint)v6 * (dpint)p5;
  c[4] = ((spint)t & mask);
  t >>= 54;
  t += (dpint)a[6] * b[6];
  t += (dpint)v6 * (dpint)p6;
  c[5] = ((spint)t & mask);
  t >>= 54;
  c[6] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
static void inline modsqr(const spint *a, spint *c) {
  udpint tot;
  udpint t = 0;
  spint p3 = 0x2488d500000000u;
  spint p4 = 0x25cfee16507a30u;
  spint p5 = 0x22c7bc80527adeu;
  spint p6 = 0xb4000a89caeau;
  spint q = ((spint)1 << 54u); // q is unsaturated radix
  spint mask = (spint)(q - (spint)1);
  tot = (udpint)a[0] * a[0];
  t = tot;
  spint v0 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  spint v1 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[2];
  tot *= 2;
  tot += (udpint)a[1] * a[1];
  t += tot;
  spint v2 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[3];
  tot += (udpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  t += (udpint)v0 * p3;
  spint v3 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[4];
  tot += (udpint)a[1] * a[3];
  tot *= 2;
  tot += (udpint)a[2] * a[2];
  t += tot;
  t += (udpint)v0 * p4;
  t += (udpint)v1 * p3;
  spint v4 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[5];
  tot += (udpint)a[1] * a[4];
  tot += (udpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (udpint)v0 * p5;
  t += (udpint)v1 * p4;
  t += (udpint)v2 * p3;
  spint v5 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[0] * a[6];
  tot += (udpint)a[1] * a[5];
  tot += (udpint)a[2] * a[4];
  tot *= 2;
  tot += (udpint)a[3] * a[3];
  t += tot;
  t += (udpint)v0 * p6;
  t += (udpint)v1 * p5;
  t += (udpint)v2 * p4;
  t += (udpint)v3 * p3;
  spint v6 = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[1] * a[6];
  tot += (udpint)a[2] * a[5];
  tot += (udpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (udpint)v1 * p6;
  t += (udpint)v2 * p5;
  t += (udpint)v3 * p4;
  t += (udpint)v4 * p3;
  c[0] = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[2] * a[6];
  tot += (udpint)a[3] * a[5];
  tot *= 2;
  tot += (udpint)a[4] * a[4];
  t += tot;
  t += (udpint)v2 * p6;
  t += (udpint)v3 * p5;
  t += (udpint)v4 * p4;
  t += (udpint)v5 * p3;
  c[1] = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[3] * a[6];
  tot += (udpint)a[4] * a[5];
  tot *= 2;
  t += tot;
  t += (udpint)v3 * p6;
  t += (udpint)v4 * p5;
  t += (udpint)v5 * p4;
  t += (udpint)v6 * p3;
  c[2] = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[4] * a[6];
  tot *= 2;
  tot += (udpint)a[5] * a[5];
  t += tot;
  t += (udpint)v4 * p6;
  t += (udpint)v5 * p5;
  t += (udpint)v6 * p4;
  c[3] = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[5] * a[6];
  tot *= 2;
  t += tot;
  t += (udpint)v5 * p6;
  t += (udpint)v6 * p5;
  c[4] = ((spint)t & mask);
  t >>= 54;
  tot = (udpint)a[6] * a[6];
  t += tot;
  t += (udpint)v6 * p6;
  c[5] = ((spint)t & mask);
  t >>= 54;
  c[6] = (spint)t;
}

// copy
static void inline modcpy(const spint *a, spint *c) {
  int i;
  for (i = 0; i < 7; i++) {
    c[i] = a[i];
  }
}

// square n times
static void modnsqr(spint *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    modsqr(a, a);
  }
}

// Calculate progenitor
static void modpro(const spint *w, spint *z) {
  spint x[7];
  spint t0[7];
  spint t1[7];
  spint t2[7];
  spint t3[7];
  spint t4[7];
  spint t5[7];
  spint t6[7];
  spint t7[7];
  spint t8[7];
  spint t9[7];
  spint t10[7];
  spint t11[7];
  spint t12[7];
  spint t13[7];
  spint t14[7];
  spint t15[7];
  spint t16[7];
  spint t17[7];
  modcpy(w, x);
  modsqr(x, z);
  modmul(x, z, t3);
  modmul(x, t3, t4);
  modmul(x, t4, t5);
  modmul(z, t5, t7);
  modmul(z, t7, t2);
  modmul(z, t2, t6);
  modmul(t3, t6, t0);
  modmul(t3, t0, t1);
  modmul(t4, t1, t14);
  modmul(t4, t14, t11);
  modmul(t0, t11, t16);
  modmul(z, t16, t9);
  modmul(z, t9, t15);
  modmul(z, t15, t0);
  modmul(t4, t0, t13);
  modmul(t4, t13, t0);
  modmul(z, t0, t10);
  modmul(t4, t10, t12);
  modmul(z, t12, t4);
  modsqr(t4, z);
  modmul(t5, z, t8);
  modmul(t0, t8, t17);
  modsqr(t8, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modmul(t8, z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modsqr(t17, t17);
  modsqr(z, z);
  modmul(t3, z, z);
  modnsqr(t17, 9);
  modmul(t14, t17, t17);
  modnsqr(t17, 9);
  modmul(t16, t17, t16);
  modnsqr(t16, 8);
  modmul(t15, t16, t15);
  modnsqr(t15, 5);
  modmul(t14, t15, t14);
  modnsqr(t14, 6);
  modmul(t1, t14, t14);
  modnsqr(t14, 7);
  modmul(t13, t14, t13);
  modnsqr(t13, 6);
  modmul(t12, t13, t12);
  modnsqr(t12, 5);
  modmul(t11, t12, t11);
  modnsqr(t11, 14);
  modmul(t9, t11, t11);
  modnsqr(t11, 8);
  modmul(t4, t11, t11);
  modnsqr(t11, 7);
  modmul(t10, t11, t10);
  modnsqr(t10, 6);
  modmul(t9, t10, t9);
  modnsqr(t9, 4);
  modmul(t7, t9, t9);
  modnsqr(t9, 9);
  modmul(t8, t9, t8);
  modnsqr(t8, 4);
  modmul(t7, t8, t7);
  modnsqr(t7, 8);
  modmul(t6, t7, t6);
  modnsqr(t6, 5);
  modmul(t5, t6, t5);
  modnsqr(t5, 11);
  modmul(t4, t5, t4);
  modnsqr(t4, 5);
  modmul(t3, t4, t3);
  modnsqr(t3, 8);
  modmul(t2, t3, t2);
  modnsqr(t2, 7);
  modmul(t1, t2, t1);
  modnsqr(t1, 9);
  modmul(t0, t1, t0);
  modnsqr(t0, 18);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, t0);
  modnsqr(t0, 16);
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
static void modinv(const spint *x, const spint *h, spint *z) {
  spint s[7];
  spint t[7];
  if (h == NULL) {
    modpro(x, t);
  } else {
    modcpy(h, t);
  }
  modcpy(x, s);
  modnsqr(t, 2);
  modmul(s, t, z);
}

// Convert m to n-residue form, n=nres(m)
static void nres(const spint *m, spint *n) {
  const spint c[7] = {0x1eb1ce624441aau, 0x8fcf6e73a8206u,  0x15642716793160u,
                      0x35056d8a688297u, 0x1982a85bf5bfa9u, 0x1bd57e579fee08u,
                      0x4e3343a49f56u};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
static void redc(const spint *n, spint *m) {
  int i;
  spint c[7];
  c[0] = 1;
  for (i = 1; i < 7; i++) {
    c[i] = 0;
  }
  modmul(n, c, m);
  (void)modfsb(m);
}

// is unity?
static int modis1(const spint *a) {
  int i;
  spint c[7];
  spint c0;
  spint d = 0;
  redc(a, c);
  for (i = 1; i < 7; i++) {
    d |= c[i];
  }
  c0 = (spint)c[0];
  return ((spint)1 & ((d - (spint)1) >> 54u) &
          (((c0 ^ (spint)1) - (spint)1) >> 54u));
}

// is zero?
static int modis0(const spint *a) {
  int i;
  spint c[7];
  spint d = 0;
  redc(a, c);
  for (i = 0; i < 7; i++) {
    d |= c[i];
  }
  return ((spint)1 & ((d - (spint)1) >> 54u));
}

// set to zero
static void modzer(spint *a) {
  int i;
  for (i = 0; i < 7; i++) {
    a[i] = 0;
  }
}

// set to one
static void modone(spint *a) {
  int i;
  a[0] = 1;
  for (i = 1; i < 7; i++) {
    a[i] = 0;
  }
  nres(a, a);
}

// set to integer
static void modint(int x, spint *a) {
  int i;
  a[0] = (spint)x;
  for (i = 1; i < 7; i++) {
    a[i] = 0;
  }
  nres(a, a);
}

// Modular multiplication by an integer, c=a*b mod 2p
// uses special method for trinomials, otherwise Barrett-Dhem reduction
static void inline modmli(const spint *a, int b, spint *c) {
  spint p3 = 0x2488d500000000u;
  spint p4 = 0x25cfee16507a30u;
  spint p5 = 0x22c7bc80527adeu;
  spint p6 = 0xb4000a89caeau;
  spint mask = ((spint)1 << 54u) - (spint)1;
  udpint t = 0;
  spint q, h, r = 0x5b05ab06ce6b27;
  t += (udpint)a[0] * (udpint)b;
  c[0] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[1] * (udpint)b;
  c[1] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[2] * (udpint)b;
  c[2] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[3] * (udpint)b;
  c[3] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[4] * (udpint)b;
  c[4] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[5] * (udpint)b;
  c[5] = (spint)t & mask;
  t = t >> 54u;
  t += (udpint)a[6] * (udpint)b;
  c[6] = (spint)t & mask;

  // Barrett-Dhem reduction
  h = (spint)(t >> 38u);
  q = (spint)(((udpint)h * (udpint)r) >> 64u);
  c[0] += q;
  t = (udpint)q * (udpint)p3;
  c[3] -= (spint)t & mask;
  c[4] -= (spint)(t >> 54u);
  t = (udpint)q * (udpint)p4;
  c[4] -= (spint)t & mask;
  c[5] -= (spint)(t >> 54u);
  t = (udpint)q * (udpint)p5;
  c[5] -= (spint)t & mask;
  c[6] -= (spint)(t >> 54u);
  c[6] = (c[6] - (q * p6)) & mask;
  (void)prop(c);
}

// Test for quadratic residue
static int modqr(const spint *h, const spint *x) {
  spint r[7];
  if (h == NULL) {
    modpro(x, r);
    modsqr(r, r);
  } else {
    modsqr(h, r);
  }
  modmul(r, x, r);
  return modis1(r) | modis0(x);
}

// conditional move g to f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__((noinline)) modcmv(int b, const spint *g,
                                             volatile spint *f) {
  int i;
  spint c0, c1, s, t, w, aux;
  static spint R = 0;
  R += 0x3cc3c33c5aa5a55au;
  w = R;
  c0 = (~b) & (w + 1);
  c1 = b + w;
  for (i = 0; i < 7; i++) {
    s = g[i];
    t = f[i];
    f[i] = aux = c0 * t + c1 * s;
    f[i] = aux - w * (t + s);
  }
}

// conditional swap g and f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__((noinline)) modcsw(int b, volatile spint *g,
                                             volatile spint *f) {
  int i;
  spint c0, c1, s, t, w, v, aux;
  static spint R = 0;
  R += 0x3cc3c33c5aa5a55au;
  w = R;
  c0 = (~b) & (w + 1);
  c1 = b + w;
  for (i = 0; i < 7; i++) {
    s = g[i];
    t = f[i];
    v = w * (t + s);
    f[i] = aux = c0 * t + c1 * s;
    f[i] = aux - v;
    g[i] = aux = c0 * s + c1 * t;
    g[i] = aux - v;
  }
}

// Modular square root, provide progenitor h if available, NULL if not
static void modsqrt(const spint *x, const spint *h, spint *r) {
  spint s[7];
  spint y[7];
  if (h == NULL) {
    modpro(x, y);
  } else {
    modcpy(h, y);
  }
  modmul(y, x, s);
  modcpy(s, r);
}

// shift left by less than a word
static void modshl(unsigned int n, spint *a) {
  int i;
  a[6] = ((a[6] << n)) | (a[5] >> (54u - n));
  for (i = 5; i > 0; i--) {
    a[i] = ((a[i] << n) & (spint)0x3fffffffffffff) | (a[i - 1] >> (54u - n));
  }
  a[0] = (a[0] << n) & (spint)0x3fffffffffffff;
}

// shift right by less than a word. Return shifted out part
static int modshr(unsigned int n, spint *a) {
  int i;
  spint r = a[0] & (((spint)1 << n) - (spint)1);
  for (i = 0; i < 6; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (54u - n)) & (spint)0x3fffffffffffff);
  }
  a[6] = a[6] >> n;
  return r;
}

// set a= 2^r
static void mod2r(unsigned int r, spint *a) {
  unsigned int n = r / 54u;
  unsigned int m = r % 54u;
  modzer(a);
  if (r >= 47 * 8)
    return;
  a[n] = 1;
  a[n] <<= m;
  nres(a, a);
}

// export to byte array
static void modexp(const spint *a, char *b) {
  int i;
  spint c[7];
  redc(a, c);
  for (i = 46; i >= 0; i--) {
    b[i] = c[0] & (spint)0xff;
    (void)modshr(8, c);
  }
}

// import from byte array
// returns 1 if in range, else 0
static int modimp(const char *b, spint *a) {
  int i, res;
  for (i = 0; i < 7; i++) {
    a[i] = 0;
  }
  for (i = 0; i < 47; i++) {
    modshl(8, a);
    a[0] += (spint)(unsigned char)b[i];
  }
  res = modfsb(a);
  nres(a, a);
  return res;
}

// determine sign
static int modsign(const spint *a) {
  spint c[7];
  redc(a, c);
  return c[0] % 2;
}

// return true if equal
static int modcmp(const spint *a, const spint *b) {
  spint c[7], d[7];
  int i, eq = 1;
  redc(a, c);
  redc(b, d);
  for (i = 0; i < 7; i++) {
    eq &= (((c[i] ^ d[i]) - 1) >> 54) & 1;
  }
  return eq;
}


/******************************************************************************
API functions calling generated code above
******************************************************************************/

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