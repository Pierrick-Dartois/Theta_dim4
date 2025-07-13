.PHONY: target gen_gf move_file clean_file
target: gen_gf move_file clean_file

gen_gf:
	python modarith/monty.py 64 0x4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

move_file:
	mv field.c src/fp_5_2_248.c

clean_file:
	rm time.c



