
CC = gcc

CFLAGS_GF = -Igf/include -Igf/test -Icommon/include -Iparams

TARGETS_FP = gf/test/test_fp.c gf/test/test_utils.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c

test_gf: test_fp_${PNAME}

test_fp_${PNAME}: ${TARGETS_FP}
	${CC} ${CFLAGS_GF} -o test_fp_${PNAME} ${TARGETS_FP}

clean_gf:
	rm test_fp_${PNAME}
