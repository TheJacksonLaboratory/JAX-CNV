
MASTER_DIR=$(shell pwd)
OBJ_DIR=$(MASTER_DIR)/obj
BIN_DIR=$(MASTER_DIR)/bin
LIB=$(MASTER_DIR)/lib

AUTOCONF = autoconf
AUTOHEADER = autoheader

CFLAGS:=-pthread
ifeq ($(mode), debug)
	CFLAGS:=$(CFLAGS) -O0 -g -DDEBUG -D_DEBUG
else
	CFLAGS:=$(CFLAGS) -mtune=native -O3 -DNDEBUG -DRELEASE
endif

CXXFLAGS:=-std=c++11 $(CFLAGS)

SUB_DIRS = $(LIB)/fastaq $(LIB)/jellyfish-2.2.6
SOURCES = main.cpp

INCLUDE= -I lib/jellyfish-2.2.6/include -I lib/fastaq/include/

all: fastaq
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o test $(SOURCES) $(INCLUDE) $(LIB)/fastaq/obj/*.o $(LIB)/jellyfish-2.2.6/lib/*.o -lz

.PHONY: all

clean:
	$(MAKE) -C $(LIB)/fastaq
	$(MAKE) -C $(LIB)/jellyfish-2.2.6
.PHONY: clean


$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

fastaq:
	@echo "- Building in fastaq"
	@$(MAKE) --no-print-directory --directory=$(LIB)/fastaq

jellyfish:
	@echo "- Building in jellyfish"
	@rm -f $(LIB)/jellyfish-2.2.6/configure
	@rm -rf $(LIB)/jellyfish-2.2.6/autom4te.cache
	@cd $(LIB)/jellyfish-2.2.6 && $(AUTOHEADER) && $(AUTOCONF) && ./configure
	$(MAKE) --no-print-directory --directory=$(LIB)/jellyfish-2.2.6
	

