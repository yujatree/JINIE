# ==============================================================================
#  Makefile for JINIE Project
#  Compiler: Intel Fortran Compiler (ifx)
# ==============================================================================

# Variables 
FC        = ifx
FFLAGS    = -O3
SRC_DIR   = src
BIN_DIR   = bin
EXAMPLE_DIR = example

PROGRAMS = D2L GRT ICD MSD NGP RDF
SRCS     = $(addsuffix .f90, $(PROGRAMS))
OBJS     = $(addsuffix .o, $(PROGRAMS))
EXES     = $(addprefix $(BIN_DIR)/, $(PROGRAMS))

# Default Target
all: dirs $(EXES)

# Directory Setup
dirs:
	@mkdir -p $(BIN_DIR)

# Module Compilation (READER)
$(SRC_DIR)/READER.mod $(SRC_DIR)/READER.o: $(SRC_DIR)/READER.f90
	$(FC) $(FFLAGS) -c $< -o $(SRC_DIR)/READER.o

# Program Compilation

$(BIN_DIR)/%: $(SRC_DIR)/%.f90 $(SRC_DIR)/READER.o
	@echo "Compiling and linking $< to $@"
	# 1. src > obj
	$(FC) $(FFLAGS) -c $< -I$(SRC_DIR) -o $(SRC_DIR)/$*.o
	# 2. linking objs
	$(FC) $(FFLAGS) $(SRC_DIR)/$*.o $(SRC_DIR)/READER.o -o $@

# Clean Targe 
clean:
	@echo "Cleaning compiled files and directories..."
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod
	rm -f $(BIN_DIR)/*
	rmdir $(BIN_DIR) 2>/dev/null || true

.PHONY: all dirs clean
