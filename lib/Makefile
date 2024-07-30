# File Names.
LIBCBRAID_SRC = utility.cpp braiding.cpp
LIBCBRAID_GROUP_SRC = artin.cpp band.cpp octahedral.cpp dihedral.cpp dual_complex.cpp standard_complex.cpp
LIBCBRAID_SRC_CLEAN = $(LIBCBRAID_SRC) $(LIBCBRAID_GROUP_SRC)
LIBCBRAID_INCLUDEDIR = ../include
LIBCBRAID_OBJ = $(LIBCBRAID_SRC:.cpp=.o)
LIBCBRAID_OBJ_CLEAN = $(LIBCBRAID_SRC_CLEAN:.cpp=.o)
LIBCBRAID_A = libcbraid.a
LIBCBRAID_CLEAN = $(LIBCBRAID_OBJ_CLEAN) $(LIBCBRAID_A)

# Compile options.
# Use make USE_XXX=1 to compile XXXxxx.cpp.
# Use make USE_ALL=1 to compile everything.

# Compile artin_braid.cpp.
ifeq ($(USE_ARTIN), 1)
	LIBCBRAID_SRC += artin.cpp
endif

# Compile band_braid.cpp.
# Band braids need CLN to be randomizable.
ifeq ($(USE_BAND), 1)
	LIBCBRAID_SRC += band.cpp
endif

ifeq ($(USE_OCTAHEDRAL), 1)
	LIBCBRAID_SRC += octahedral.cpp
endif
ifeq ($(USE_DIHEDRAL), 1)
	LIBCBRAID_SRC += dihedral.cpp
endif
ifeq ($(USE_DUAL_COMPLEX), 1)
	LIBCBRAID_SRC += dual_complex.cpp
endif
ifeq ($(USE_STANDARD_COMPLEX), 1)
	LIBCBRAID_SRC += standard_complex.cpp
endif
ifeq ($(USE_ALL),1)
	LIBCBRAID_SRC += $(LIBCBRAID_GROUP_SRC)
endif
# Use make USE_CLN=1 to compile with the CLN library.
ifeq ($(USE_CLN),1)
	CPPFLAGS_CLN = -DUSE_CLN
	LIBFLAGS_CLN = -lcln
endif
CPPFLAGS = -Wall -O3 -I$(LIBCBRAID_INCLUDEDIR) $(CPPFLAGS_CLN) -DBRAIDING_CLASS=$(BRAIDING_CLASS) $(CPPFLAGS_PAR)
LIBFLAGS = -lcbraid -lm $(LIBFLAGS_CLN)
ARFLAGS = rv
CXX = g++

MAKEDEPCMD = $(CXX)
MAKEDEPFLAGS = -M
DEPFILE = .depends

# Default target.
all: $(LIBCBRAID_A)

# Default compile rules.
.cpp.o:
	$(CXX) $(CPPFLAGS) -c $<

# Build libcbraid.
$(LIBCBRAID_A): $(LIBCBRAID_OBJ)
	$(AR) $(ARFLAGS) $@ $?

$(LIBCBRAID_OBJ): $(LIBCBRAID_SRC)

# Cleanup.
clean:
	rm -f $(LIBCBRAID_CLEAN)
	rm -f $(DEPFILE)

# Make dependencies.
depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

# Make dependencies file.
$(DEPFILE):
	$(MAKEDEPCMD) $(MAKEDEPFLAGS) $(LIBCBRAID_SRC) \
		-I$(LIBCBRAID_INCLUDEDIR) > $(DEPFILE)

# Read in dependencies file created by 'make depend'.
sinclude $(DEPFILE)
