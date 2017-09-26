
MASTER_DIR=$(shell pwd)
OBJ_DIR=$(MASTER_DIR)/obj
BIN_DIR=$(MASTER_DIR)/bin
LIB=$(MASTER_DIR)/lib

AUTOCONF = autoconf
AUTOHEADER = autoheader

CFLAGS:=-pthread
ifeq ($(mode), debug)
	CFLAGS:=$(CFLAGS) -O0 -g
else
	CFLAGS:=$(CFLAGS) -mtune=native -O3 -Wall
endif

CXXFLAGS:=-std=c++11 $(CFLAGS)

SUB_DIRS = $(LIB)/fastaq $(LIB)/jellyfish-2.2.6
SOURCES = main.cpp

PROGRAM=$(BIN_DIR)/GetKmerCount
JELLYFISH=$(LIB)/jellyfish-2.2.6/bin/jellyfish

INCLUDE= -I lib/jellyfish-2.2.6/include -I lib/fastaq/include/ -I lib/htslib/

all: $(PROGRAM)
.PHONY: all

$(PROGRAM): fastaq $(JELLYFISH) $(SOURCES)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(INCLUDE) $(LIB)/fastaq/obj/*.o $(LIB)/jellyfish-2.2.6/lib/*.o -lz

.PHONY: all

clean:
	$(MAKE) clean -C $(LIB)/fastaq
	@rm -rf $(LIB)/jellyfish-2.2.6
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
.PHONY: clean


$(OBJ_DIR):
	@mkdir -p $@
$(BIN_DIR):
	@mkdir -p $@

fastaq:
	@echo "- Building in fastaq"
	@$(MAKE) --no-print-directory --directory=$(LIB)/fastaq

$(JELLYFISH): $(BIN_DIR)
	@echo "- Building in jellyfish"
	@cd $(LIB) && tar -zxvf $(LIB)/jellyfish-2.2.6.tar.gz
	@cd $(LIB)/jellyfish-2.2.6 && ./configure --prefix=$(LIB)/jellyfish-2.2.6
	$(MAKE) --no-print-directory --directory=$(LIB)/jellyfish-2.2.6
	@cp $@ $(BIN_DIR)

