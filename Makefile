
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

INCLUDE= -I lib/jellyfish-2.2.6/include -I lib/fastaq/include/

all: $(PROGRAM)
.PHONY: all

$(PROGRAM): fastaq jellyfish
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(INCLUDE) $(LIB)/fastaq/obj/*.o $(LIB)/jellyfish-2.2.6/lib/*.o -lz

.PHONY: all

clean:
	$(MAKE) clean -C $(LIB)/fastaq
	$(MAKE) clean -C $(LIB)/jellyfish-2.2.6
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
.PHONY: clean


$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

fastaq:
	@echo "- Building in fastaq"
	@$(MAKE) --no-print-directory --directory=$(LIB)/fastaq

jellyfish:
	@echo "- Building in jellyfish"
	@cd $(LIB)/jellyfish-2.2.6 && ./configure --prefix=$(LIB)/jellyfish-2.2.6
	$(MAKE) --no-print-directory --directory=$(LIB)/jellyfish-2.2.6
	

