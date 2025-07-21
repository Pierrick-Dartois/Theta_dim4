#include <intbig.h> // Where ibz_t is defined
#define RADIX 64
extern const uint64_t NWORDS_FIELD;
extern const uint64_t NWORDS_ORDER;
extern const uint64_t TORSION_EVEN_POWER;//e_2
extern const ibz_t CHARACTERISTIC;
extern const ibz_t TORSION_EVEN;// 2^e_2
extern const ibz_t TORSION_ODD;// c*l^e_l
#if defined(TEST_KANI)
extern const uint64_t L;// l
extern const uint64_t TORSION_L_POWER;// e_l
extern const ibz_t TORSION_L;// l^e_l
extern const ibz_t TORSION_COMPLEMENT_L;// c*2^e_2
// a1**2+a2**2+l**fl=2**f2
extern const uint64_t CONST_FL;
extern const uint64_t CONST_F2;
extern const ibz_t CONST_A1;
extern const ibz_t CONST_A2;

#elif defined(SIDH_ATTACK)
extern const uint64_t L;// l
extern const uint64_t TORSION_L_POWER;// e_l
extern const ibz_t TORSION_L;// l^e_l
extern const ibz_t TORSION_COMPLEMENT_L;// c*2^e_2
// a1**2+a2**2+l**fl=2**f2
extern const uint64_t CONST_FL;
extern const uint64_t CONST_F2;
extern const ibz_t CONST_A1;
extern const ibz_t CONST_A2;

#endif
