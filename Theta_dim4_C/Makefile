CODE_DIR = fiat-crypto

.PHONY: compile_fiat_crypto

FIAT_CRYPTO_BIN = fiat-crypto/src/ExtractionOcaml/bedrock2_fiat_crypto

ifneq ("$(wildcard $(FIAT_CRYPTO_BIN))","")
    MSG_FIAT_CRYPTO = Fiat-crypto is already compiled.
    FIAT_CRYPTO_SCRIPT = 
else
    MSG_FIAT_CRYPTO = Compiling Fiat-crypto...
    FIAT_CRYPTO_SCRIPT = $(MAKE) -C $(CODE_DIR) standalone-ocaml
endif

compile_fiat_crypto:
	@echo $(MSG_FIAT_CRYPTO)
	$(FIAT_CRYPTO_SCRIPT)

#all:
	#@echo $(MSG_FIAT_CRYPTO)
	#@echo $(FIAT_CRYPTO_SCRIPT)