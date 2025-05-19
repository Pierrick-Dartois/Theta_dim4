from sage.all import *
proof.all(False) 

import linecache
#from Cryptodome.Hash import SHAKE256

from utilities.basis_from_hints import torsion_basis_2f_from_hint
from utilities.discrete_log import weil_pairing_pari, discrete_log_pari
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only

public_params = {1:{'p':5*2**248-1,'c':5,'e':248},3:{'p':65*2**376-1,'c':65,'e':376},5:{'p':27*2**500-1,'c':27,'e':500}}


class SQIsignHD:
	def __init__(self,lvl):
		self.p = public_params[lvl]['p']
		self.c = public_params[lvl]['c']
		self.e = public_params[lvl]['e']
		self.Fp2 = GF(self.p**2,'i',modulus=[1,0,1],proof=False)
		self.i = self.Fp2.gen()

		self.NQR_TABLE = self.read_gf_table("API/SQIsignHD/Data/NQR_TABLE_lvl"+str(lvl)+".txt")
		self.Z_NQR_TABLE = self.read_gf_table("API/SQIsignHD/Data/Z_NQR_TABLE_lvl"+str(lvl)+".txt")

		for x in self.NQR_TABLE:
			assert(x.is_square() == False)
		for x in self.Z_NQR_TABLE:
			assert(x.is_square() and (x-1).is_square() == False)

		self.lvl = lvl
		self.lamb = 96+32*lvl
		self.n_bytes = 2*self.lamb//8
		self.f = self.lamb + ceil(log(2*self.lamb)/log(2))
		self.r = ceil(self.f/2)+2


	def read_gf_table(self,file):
		L=[]
		for j in range(20):
			line = linecache.getline(file,j+1)
			M = line.split('+')
			x = int(M[0][:-1],16)+self.i*int(M[1][3:],16)
			L.append(x)
		return L

	def read_public_key(self,file,number):
		# Starting from 0
		n_line=number*3+1
	
		# A_pk
		line = linecache.getline(file,n_line)
		L=line.split('=')
		L=L[1][1:-1].split('+')
		A_pk=int(L[0][:-1],16)+self.i*int(L[1][3:],16)

		# hints
		line = linecache.getline(file,n_line+1)
		L=line.split('=')
		hint_pk_P = int(L[1][1:-1])

		line = linecache.getline(file,n_line+2)
		L=line.split('=')
		hint_pk_Q = int(L[1][1:-1])

		return A_pk,hint_pk_P,hint_pk_Q

	def read_signature(self,file,number):
		# number starting from 0
		n_line=number*11+1

		# A_com
		line = linecache.getline(file,n_line)
		L=line.split('=')
		L=L[1][1:-1].split('+')
		A_com=int(L[0][:-1],16)+self.i*int(L[1][3:],16)

		# a
		line = linecache.getline(file,n_line+1)
		L=line.split('=')
		a = int('0x'+L[1][1:-1],16)

		# b
		line = linecache.getline(file,n_line+2)
		L=line.split('=')
		b = int('0x'+L[1][1:-1],16)

		# c_or_d
		line = linecache.getline(file,n_line+3)
		L=line.split('=')
		c_or_d = int('0x'+L[1][1:-1],16)

		# q
		line = linecache.getline(file,n_line+4)
		L=line.split('=')
		q = int('0x'+L[1][1:-1],16)

		# h_com_P
		line = linecache.getline(file,n_line+5)
		L=line.split('=')
		hint_com_P = int(L[1][1:-1])

		# h_com_Q
		line = linecache.getline(file,n_line+6)
		L=line.split('=')
		hint_com_Q = int(L[1][1:-1])

		# h_chal_P
		line = linecache.getline(file,n_line+7)
		L=line.split('=')
		hint_chal_P = int(L[1][1:-1])

		# h_chal_Q
		line = linecache.getline(file,n_line+8)
		L=line.split('=')
		hint_chal_Q = int(L[1][1:-1])

		# vec_chal
		vec_chal = [0,0]
		line = linecache.getline(file,n_line+9)
		L=line.split('=')
		vec_chal[0] = int('0x'+L[1][1:-1],16)
		line = linecache.getline(file,n_line+10)
		L=line.split('=')
		vec_chal[1] = int('0x'+L[1][1:-1],16)


		return A_com, a, b, c_or_d, q, hint_com_P, hint_com_Q, hint_chal_P, hint_chal_Q, vec_chal

class SQIsignHD_verif:
	def __init__(self,pp,number):
		self.params = pp

		file = "API/SQIsignHD/Data/Public_keys_lvl"+str(pp.lvl)+".txt"
		self.A_pk, self.h_pk_P, self.h_pk_Q = pp.read_public_key(file,number)

		file = "API/SQIsignHD/Data/Signatures_lvl"+str(pp.lvl)+".txt"
		self.A_com, self.a, self.b, self.c_or_d, self.q, self.h_com_P, self.h_com_Q,\
		self.h_chal_P, self.h_chal_Q, self.vec_chal = pp.read_signature(file,number)


	def recover_pk_and_com(self):
		self.E_pk = EllipticCurve([0,self.A_pk,0,1,0])
		self.P_pk, self.Q_pk = torsion_basis_2f_from_hint(self.E_pk,self.h_pk_P,self.h_pk_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

		self.E_com = EllipticCurve([0,self.A_com,0,1,0])
		self.P_com, self.Q_com = torsion_basis_2f_from_hint(self.E_com,self.h_com_P,self.h_com_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

	def recover_chal(self):
		rescale = ZZ(2**(self.params.e-self.params.lamb))
		deg = ZZ(2**(self.params.lamb))

		B_pk_lamb = (rescale*self.P_pk, rescale*self.Q_pk)
		phi_chal, self.E_chal = isogeny_from_scalar_x_only(self.E_pk, deg, self.vec_chal, B_pk_lamb)

		self.P_chal, self.Q_chal = torsion_basis_2f_from_hint(self.E_chal,self.h_chal_P,self.h_chal_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

	def image_response(self):
		rescale = ZZ(22**(self.params.e-self.params.r))
		R_com = rescale*self.P_com
		S_com = rescale*self.Q_com
		R_chal = rescale*self.P_chal
		S_chal = rescale*self.Q_chal
		order = ZZ(2**self.params.r)

		w_com = weil_pairing_pari(R_com, S_com, order)
		w_chal = weil_pairing_pari(R_chal, S_chal, order)

		k = discrete_log_pari(w_com, w_chal, order)

		if self.a%2==1:
			self.c = self.c_or_d
			self.d = (inverse_mod(ZZ(self.a),order)*(k*self.q+self.b*self.c))%order
		else:
			self.d = self.c_or_d
			self.c = (inverse_mod(ZZ(self.b),order)*(self.a*self.d-k*self.q))%order
		self.phi_rsp_R_com = self.a*R_chal + self.b*S_chal
		self.phi_rsp_S_com = self.c*R_chal + self.d*S_chal



if __name__=="__main__":
	pp = SQIsignHD(1)
		
	for i in range(100):
		print(i)
		verif = SQIsignHD_verif(pp,i)
		verif.recover_pk_and_com()
		verif.recover_chal()
		verif.image_response()
		
