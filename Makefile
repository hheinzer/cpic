# CPIC - Particle in Cell Method, written in C++
# See LICENSE file for copyright and license details.
TARGET = libcpic
TEST = box_dsmc_nanbu_ei

# config {on, off}
DEBUGGING = off
NDEBUG = off
OPENMP = off
PROFILING = off

# programs
CC = g++
AR = ar rcs

# source file ending
C = cpp

.PHONY: all clean run debug memdebug

# libs and incs
LIBS =
INCS = $(shell find src -type d -exec echo -I{} \;)

# flags
FLAGS = -Wall -Wextra -pedantic -pipe -ggdb3
ifeq ($(DEBUGGING), on)
  FLAGS += -O0
else
  FLAGS += -O3 -march=native -flto
  ifeq ($(NDEBUG), on)
    FLAGS += -DNDEBUG
  endif
endif
ifeq ($(OPENMP), on)
  FLAGS += -fopenmp
else
  FLAGS += -Wno-unknown-pragmas
endif
ifeq ($(PROFILING), on)
  FLAGS += -pg
endif

# sources, objects, and target
SRC = $(shell find src -type f -name *.$(C))
OBJ = $(patsubst src/%.$(C), obj/%.o, $(SRC))
DEP = $(patsubst %.o, %.d, $(OBJ))
BIN = bin/$(TEST)
LIB = lib/$(TARGET)

# default
all: $(LIB) $(BIN)
lib: $(LIB)

# dependencies
DEPFLAGS = -MMD -MF $(@:.o=.d)
-include $(DEP)

# recipes
$(BIN): bin/%: test/%.$(C) $(LIB)
	@mkdir -p $(@D)
	$(CC) $(FLAGS) $(INCS) $< $(LIB) -o $@ $(LIBS)

$(LIB): $(OBJ)
	@mkdir -p $(@D)
	$(AR) $@ $^

$(OBJ): obj/%.o: src/%.$(C) Makefile
	@mkdir -p $(@D)
	$(CC) $(FLAGS) $(INCS) -c $< -o $@ $(DEPFLAGS)

clean:
	rm -rf obj bin lib

run:
	@./$(BIN)

debug:
	@gdb ./$(BIN)

memdebug:
	@valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(BIN)
