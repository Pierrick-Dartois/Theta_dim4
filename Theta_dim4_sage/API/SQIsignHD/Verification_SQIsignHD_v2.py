from sage.all import *
proof.all(False) 

import linecache

from utilities.basis_from_hints import torsion_basis_2f_from_hint

public_params = {1:{'p':5*2**248-1,'c':5,'e':248},3:{'p':65*2**376-1,'c':65,'e':376},5:{'p':27*2**500-1,'c':27,'e':500}}

class SQIsignHD:
	def __init__(self,lvl):
		self.p=public_params[lvl]['p']
		self.c=public_params[lvl]['c']
		self.e=public_params[lvl]['e']
		self.Fp2=GF(self.p**2,'i',modulus=[1,0,1],proof=False)
		self.i=self.Fp2.gen()

		self.NQR_TABLE = self.read_gf_table("API/SQIsignHD/Data/NQR_TABLE_lvl"+str(lvl)+".txt")
		self.Z_NQR_TABLE = self.read_gf_table("API/SQIsignHD/Data/Z_NQR_TABLE_lvl"+str(lvl)+".txt")

		for x in self.NQR_TABLE:
			assert(x.is_square()==False)
		for x in self.Z_NQR_TABLE:
			assert(x.is_square() and (x-1).is_square()==False)


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
		# Starting from 0
		n_line=number*9+1

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

		return A_com, a, b, c_or_d, q, hint_com_P, hint_com_Q, hint_chal_P, hint_chal_Q

if __name__=="__main__":
	inst = SQIsignHD(1)
	for i in range(100):
		A,hP,hQ = inst.read_public_key("API/SQIsignHD/Data/Public_keys_lvl1.txt",i)
		E = EllipticCurve([0,A,0,1,0])
		P,Q = torsion_basis_2f_from_hint(E,hP,hQ,inst.NQR_TABLE,inst.Z_NQR_TABLE)
		A_com, a, b, c_or_d, q, hint_com_P, hint_com_Q, hint_chal_P, hint_chal_Q = inst.read_signature("API/SQIsignHD/Data/Signatures_lvl1.txt",i)
		E_com = EllipticCurve([0,A_com,0,1,0])
		P_com,Q_com = torsion_basis_2f_from_hint(E_com,hint_com_P,hint_com_Q,inst.NQR_TABLE,inst.Z_NQR_TABLE)

