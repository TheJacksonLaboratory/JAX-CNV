
MASTER_DIR=$(shell pwd)
OBJ_DIR=$(MASTER_DIR)/obj
BIN_DIR=$(MASTER_DIR)/bin
LIB=$(MASTER_DIR)/lib

AUTOCONF = autoconf
AUTOHEADER = autoheader

CFLAGS:=
ifeq ($(mode), debug)
	CFLAGS:=$(CFLAGS) -O0 -g -DDEBUG -D_DEBUG
else
	CFLAGS:=$(CFLAGS) -mtune=native -O3 -DNDEBUG -DRELEASE
endif

CXXFLAGS:=-std=c++11 $(CFLAGS)

SUB_DIRS = $(LIB)/fastaq
SOURCES = main.cpp

all: AUX
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o test $(SOURCES) $(LIB)/fastaq/obj/*.o -lz

.PHONY: all

clean:
.PHONY: clean


$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

AUX: $(OBJ_DIR)
	@for dir in $(SUB_DIRS); do \
		$(MAKE) --no-print-directory --directory=$$dir; \
	done

