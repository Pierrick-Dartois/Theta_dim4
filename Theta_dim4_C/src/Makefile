
CC = gcc

V_DEFS = -DP_${PNAME}

CFLAGS_GF = -Igf/include -Igf/test -Icommon/include -Iparams -Imini-gmp/include

TARGETS_FP = gf/test/test_fp.c gf/test/test_utils.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c

test_gf: test_fp_${PNAME}

test_fp_${PNAME}: ${TARGETS_FP}
	${CC} ${V_DEFS} ${CFLAGS_GF} -o test_fp_${PNAME} ${TARGETS_FP}

clean_gf:
	rm test_fp_${PNAME}
