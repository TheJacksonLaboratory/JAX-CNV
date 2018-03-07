
MASTER_DIR=$(shell pwd)
OBJ_DIR=$(MASTER_DIR)/obj
BIN_DIR=$(MASTER_DIR)/bin
LIB=$(MASTER_DIR)/lib

AUTOCONF = autoconf
AUTOHEADER = autoheader

CFLAGS:=-pthread
ifeq ($(mode), debug)
	CFLAGS:=$(CFLAGS) -O0 -g -DDEBUG
else
	CFLAGS:=$(CFLAGS) -O3
endif

CXXFLAGS:=-std=c++11 $(CFLAGS)
export $(CFLAGS)
export $(CXXFLAGS)

SOURCES = main.cpp src/GrabJellyfishKmer.cpp src/GetCnvSignal.cpp src/GenerateKmer.cpp src/CallHmm.cpp src/EstimateCoverage.cpp

PROGRAM=$(BIN_DIR)/clia_cnv

UMDHMM_SRC= $(LIB)/umdhmm-v1.02/backward.c \
		$(LIB)/umdhmm-v1.02/baum.c \
		$(LIB)/umdhmm-v1.02/forward.c \
		$(LIB)/umdhmm-v1.02/hmmrand.c \
		$(LIB)/umdhmm-v1.02/hmmutils.c \
		$(LIB)/umdhmm-v1.02/nrutil.c \
		$(LIB)/umdhmm-v1.02/sequence.c \
		$(LIB)/umdhmm-v1.02/viterbi.c

INCLUDE = -I lib/jellyfish-2.2.6/include -I lib/fastaq/include/ -I lib/ -I lib/htslib/ -I include/
LIBRARY = -lz -lcurl -lbz2 $(LIB)/fastaq/obj/*.o $(LIB)/jellyfish-2.2.6/lib/*.o \
		$(patsubst %.c, %.o, $(UMDHMM_SRC) )

JELLYFISH=$(LIB)/jellyfish-2.2.6/bin/jellyfish
HTS_LIB=$(LIB)/htslib/libhts.a

all: $(PROGRAM)
.PHONY: all

$(PROGRAM): fastaq umdhmm $(JELLYFISH) $(HTS_LIB) $(SOURCES)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(UMDHMM_SRC) $(INCLUDE) $(HTS_LIB) $(LIBRARY)

.PHONY: all

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	$(MAKE) clean -C $(LIB)/fastaq
	@rm -rf $(LIB)/jellyfish-2.2.6
	$(MAKE) clean -C $(LIB)/umdhmm-v1.02
.PHONY: clean


$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

fastaq:
	@echo "- Building in fastaq"
	@$(MAKE) --no-print-directory --directory=$(LIB)/fastaq

umdhmm:
	@echo "- Building in umdhmm"
	@$(MAKE) --no-print-directory --directory=$(LIB)/umdhmm-v1.02

$(JELLYFISH):
	@echo "- Building in jellyfish"
	@cd $(LIB) && tar -zxvf $(LIB)/jellyfish-2.2.6.tar.gz
	@cd $(LIB)/jellyfish-2.2.6 && ./configure --prefix=$(LIB)/jellyfish-2.2.6
	$(MAKE) --no-print-directory --directory=$(LIB)/jellyfish-2.2.6
	@mkdir -p $(BIN_DIR)
	@cp $(JELLYFISH) $(BIN_DIR)

$(HTS_LIB):
	@echo "- Building in htslib"
	@rm -f $(LIB)/htslib/configure
	@rm -rf $(LIB)/htslib/autom4te.cache
	@cd $(LIB)/htslib && $(AUTOHEADER) && $(AUTOCONF) && ./configure --disable-lzma --disable-lcurl
	$(MAKE) --no-print-directory -C $(LIB)/htslib
