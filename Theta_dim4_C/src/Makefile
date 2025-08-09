
CC = gcc

V_DEFS = -DP_${PNAME} -DNDEBUG

CFLAGS_GF = -Igf/include -Igf/test -Icommon/include -Iparams -Imini-gmp/include

TARGETS_TEST_FP = gf/test/test_fp.c gf/test/test_utils.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c
TARGETS_BENCH_FP = gf/test/bench_fp.c gf/test/test_utils.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c
TARGETS_TEST_FP2 = gf/test/test_fp2.c gf/test/test_utils.c gf/fp2.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c
TARGETS_BENCH_FP2 = gf/test/bench_fp2.c gf/test/test_utils.c gf/fp2.c gf/fp.c gf/fp/fp_${PNAME}.c params/fp_params/constants_${PNAME}.c common/randombytes_system.c

test_and_bench_gf: test_gf bench_gf

test_gf: test_fp_${PNAME} test_fp2_${PNAME}

test_fp_${PNAME}: ${TARGETS_TEST_FP}
	${CC} ${V_DEFS} ${CFLAGS_GF} -o test_fp_${PNAME} ${TARGETS_TEST_FP}

test_fp2_${PNAME}: ${TARGETS_TEST_FP2}
	${CC} ${V_DEFS} ${CFLAGS_GF} -o test_fp2_${PNAME} ${TARGETS_TEST_FP2}

bench_gf: bench_fp_${PNAME} bench_fp2_${PNAME}

bench_fp_${PNAME}: ${TARGETS_BENCH_FP}
	${CC} ${V_DEFS} ${CFLAGS_GF} -o bench_fp_${PNAME} ${TARGETS_BENCH_FP}

bench_fp2_${PNAME}: ${TARGETS_BENCH_FP2}
	${CC} ${V_DEFS} ${CFLAGS_GF} -o bench_fp2_${PNAME} ${TARGETS_BENCH_FP2}

clean_gf:
	rm test_fp_${PNAME} test_fp2_${PNAME} bench_fp_${PNAME} bench_fp2_${PNAME}
