from sage.all import *

import os
import argparse


L_primes=prime_range(1000)

def lazy_factor(N):
	r"""
	We return the factorization of N if it is easy 
	to factor and can be decomposed as a sum of two squares 
	and False otherwise.
	"""

	M=N
	L_fact=[]
	for p in L_primes:
		e_p=0
		while M%p==0:
			M=M//p
			e_p+=1
		if e_p>0:
			if p==2 or p%4==1 or (p%4==3 and e_p%2==0):
				L_fact.append((p,e_p))
			else:
				return False
	if is_prime(M) and M%4==1:
		L_fact.append((M,1))
		return L_fact
	else:
		return False

def lazy_cornacchia(N):
	r"""
	Returns integers a1, a2 such that N=a1**2+a2**2 if they exist and are easy to find.
	Returns False otherwise.
	"""
	Q0=BinaryQF([1,0,1])

	L_fact=lazy_factor(N)
	if L_fact:
		L_a1a2=[]
		sq_factor=1
		for fact in L_fact:
			p,e_p=fact
			if e_p%2==0:
				sq_factor *= p**(e_p//2)
			else:
				a1,a2=Q0.solve_integer(p, algorithm='cornacchia')
				L_a1a2.append((a1,a2,e_p))

		gauss_integer=1
		for elt in L_a1a2:
			a1,a2,e_p=elt
			gauss_integer*=(a1+a2*I)**e_p
		return sq_factor*gauss_integer[0],sq_factor*gauss_integer[1]
	else:
		return False

def find_embedding_params(e2_max,l,el_max):
	r"""
	Given exponents e2_max and el_max of the accessible 2**e2_max and 
	l**el_max-torsion respectively, this function finds integers e2, el, a1, a2 
	such that 2**e2=l**el+a1**2+a2**2, where el is the biggest possible.
	"""

	e2 = e2_max
	el = el_max

	for el in range(el_max,0,-1):
		for e2 in range(e2_max):
			N = 2**e2 - l**el
			if N>0:
				a1a2 = lazy_cornacchia(N)
				if a1a2:
					a1,a2 = a1a2
					return e2,el,a1,a2

	return False

def int_to_montgemery_fp_const(x,p,Nlimbs,Radix):
	R = 2**(Nlimbs*Radix)
	el = (x*R)%p
	vs = [(int(el) >> Radix*i) % 2**Radix for i in range(Nlimbs)]
	return '{' + ', '.join(map(hex, vs)) + '}'


def write_field_file(p):
	d_word_params = {}
	file = open("field.c","a+")
	file.seek(0)

	L_file = file.readlines()
	for i in range(14,19):
		L_line = L_file[i][:-1].split(" ")
		d_word_params[L_line[1]] = int(L_line[2])

	Nlimbs = d_word_params['Nlimbs']
	Radix = d_word_params['Radix']

	file.write("\n\n/******************************************************************************")
	file.write("\nAPI functions calling generated code above")
	file.write("\n******************************************************************************/")

	file.write("\n\nconst digit_t ZERO[NWORDS_FIELD] = "+int_to_montgemery_fp_const(0,p,Nlimbs,Radix)+";")
	file.write("\n\nconst digit_t ONE[NWORDS_FIELD] = "+int_to_montgemery_fp_const(1,p,Nlimbs,Radix)+";")
	file.write("\n\nstatic const digit_t TWO_INV[NWORDS_FIELD] = "+int_to_montgemery_fp_const(inverse_mod(2,p),p,Nlimbs,Radix)+";")
	file.write("\n\nstatic const digit_t THREE_INV[NWORDS_FIELD] = "+int_to_montgemery_fp_const(inverse_mod(3,p),p,Nlimbs,Radix)+";")


	file.write("\n\nvoid")
	file.write("\nfp_set_small(fp_t *x, const digit_t val)")
	file.write("\n{")
	file.write("\n    modint((int)val, *x);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)")
	file.write("\n{")
	file.write("\n    modmli(*a, (int)val, *x);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_set_zero(fp_t *x)")
	file.write("\n{")
	file.write("\n    modzer(*x);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_set_one(fp_t *x)")
	file.write("\n{")
	file.write("\n    modone(*x);")
	file.write("\n}")

	file.write("\n\nuint32_t")
	file.write("\nfp_is_equal(const fp_t *a, const fp_t *b)")
	file.write("\n{")
	file.write("\n    return -(uint32_t)modcmp(*a, *b);")
	file.write("\n}")

	file.write("\n\nuint32_t")
	file.write("\nfp_is_zero(const fp_t *a)")
	file.write("\n{")
	file.write("\n    return -(uint32_t)modis0(*a);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_copy(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modcpy(*a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_cswap(fp_t *a, fp_t *b, uint32_t ctl)")
	file.write("\n{")
	file.write("\n    modcsw((int)(ctl & 0x1), *a, *b);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_add(fp_t *out, const fp_t *a, const fp_t *b)")
	file.write("\n{")
	file.write("\n    modadd(*a, *b, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_sub(fp_t *out, const fp_t *a, const fp_t *b)")
	file.write("\n{")
	file.write("\n    modsub(*a, *b, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_neg(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modneg(*a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_sqr(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modsqr(*a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_mul(fp_t *out, const fp_t *a, const fp_t *b)")
	file.write("\n{")
	file.write("\n    modmul(*a, *b, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_inv(fp_t *x)")
	file.write("\n{")
	file.write("\n    modinv(*x, NULL, *x);")
	file.write("\n}")

	file.write("\n\nuint32_t")
	file.write("\nfp_is_square(const fp_t *a)")
	file.write("\n{")
	file.write("\n    return -(uint32_t)modqr(NULL, *a);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_sqrt(fp_t *a)")
	file.write("\n{")
	file.write("\n    modsqrt(*a, NULL, *a);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_half(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modmul(TWO_INV, *a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_exp3div4(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modpro(*a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("\nfp_div3(fp_t *out, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modmul(THREE_INV, *a, *out);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("fp_encode(void *dst, const fp_t *a)")
	file.write("\n{")
	file.write("\n    modexp(a, dst);")
	file.write("\n}")

	file.write("\n\nvoid")
	file.write("fp_decode_reduce(fp_t *d, const void *src, size_t len)")
	file.write("\n{")
	file.write("\n    modimp(src,d);")
	file.write("\n}")






if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-p","--prime")
	parser.add_argument("-n","--name")
	parser.add_argument("-t","--test",action="store_true")

	args = parser.parse_args()

	p = eval(args.prime)
	L_factors = list(factor(p+1))
	power_2 = L_factors[0][1]

	if args.test:
		# Find l, power_l such that p = c * 2**power_2 * l**power_l - 1
		i = 1
		power_l = L_factors[1][1]
		for j in range(1,len(L_factors)):
			if L_factors[j][1] > power_l:
				i = j

		l = L_factors[i][0]
		power_l = L_factors[i][1]

		# Find f_2, f_l, a1, a2 such that f_l <= power_l, f_2 <= power_2-2 and a_1**2 + a2**2 + l**f_l = 2**f_2
		f_2, f_l, a1, a2 = find_embedding_params(power_2,l,power_l)

	os.system("python ../modarith/monty.py 64 {}".format(p))

	write_field_file(p)

	os.system("mv field.c ../src/gf/fp/fp_"+args.name+".c")
	os.system("rm time.c")