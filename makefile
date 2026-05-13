# ==============================================================================
#  Makefile for JINIE Project
#  Compiler: Intel Fortran Compiler (ifx)
# ==============================================================================

FC        = ifx
FFLAGS    = -O3
LD_FLAGS  =
SRC_DIR   = src
BIN_DIR   = bin

# Program groups----------------------------------------------------------------

PROGRAMS_READER = D2L AGD CLS CLF RDF R2X R2N MSD NGP ICD GSRT RGRT HOP h_vis CLEAR_TRAJ HELP_EFF SIEVE_TRAJ HELP_assign# #RNGP

PROGRAMS_ENSEMBLE = CRDF CDDF_Li CDDF_P CDDF_S CST

PROGRAMS = $(PROGRAMS_READER) $(PROGRAMS_ENSEMBLE)

EXES = $(addprefix $(BIN_DIR)/,$(PROGRAMS))

# Default target----------------------------------------------------------------

all: dirs $(EXES)

dirs:
	@mkdir -p $(BIN_DIR)

# READER based program----------------------------------------------------------

$(BIN_DIR)/%: $(SRC_DIR)/%.f90
	@echo "Building $@"
	@if echo "$(PROGRAMS_READER)" | grep -w "$*" > /dev/null; then \
		$(FC) $(FFLAGS) -module $(BIN_DIR) -c $(SRC_DIR)/READER.f90 -o $(BIN_DIR)/READER.o; \
		$(FC) $(FFLAGS) -module $(BIN_DIR) $(BIN_DIR)/READER.o $< -o $@; \
		rm -f $(BIN_DIR)/READER.o; \
	elif echo "$(PROGRAMS_ENSEMBLE)" | grep -w "$*" > /dev/null; then \
		$(FC) $(FFLAGS) -module $(BIN_DIR) -c $(SRC_DIR)/ENSEMBLE_READER.f90 -o $(BIN_DIR)/ENSEMBLE_READER.o; \
		$(FC) $(FFLAGS) -module $(BIN_DIR) $(BIN_DIR)/ENSEMBLE_READER.o $< -o $@; \
		rm -f $(BIN_DIR)/ENSEMBLE_READER.o; \
	else \
		echo "Error: Unknown program $*"; \
		exit 1; \
	fi

# ------------------------------------------------------------------------------
# Clean
# ------------------------------------------------------------------------------

clean:
	@echo "Cleaning..."
	rm -rf $(BIN_DIR)

.PHONY: all clean dirs

