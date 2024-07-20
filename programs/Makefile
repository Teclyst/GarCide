# File names and directories
SOURCES       = test.cpp speedtest.cpp braiding_main.cpp compare_test.cpp
PROGRAMS      = $(SOURCES:.cpp=)
OBJS	      = $(SOURCES:.cpp=.o)
CBRAID_LIBDIR = ../lib
CBRAID_INCDIR = ../include
CBRAID_LIB    = $(CBRAID_LIBDIR)/libcbraid.a
DEPFILE       = .depends

# Commands
CPP          = g++
MAKEDEPCMD   = $(CPP)

# Compile options.
# Use make USE_CLN=1 to compile with the CLN library.

# Use parallelism by default.
USE_PAR = 1

ifeq ($(USE_ALL),1)
	CPPFLAGS_CLN = -DUSE_CLN
	LIBFLAGS_CLN = -lcln
endif
ifeq ($(USE_BAND),1)
	CPPFLAGS_CLN = -DUSE_CLN
	LIBFLAGS_CLN = -lcln
endif
ifeq ($(USE_CLN),1)
	CPPFLAGS_CLN = -DUSE_CLN
	LIBFLAGS_CLN = -lcln
endif
# Compile with parallel algorithms.
# To disable parallelism, compile with USE_PAR=0 (or anything else than 1).
ifeq ($(USE_PAR),1)
	CPPFLAGS_PAR = -DUSE_PAR
endif
ifeq ($(RANDOMIZE_ON_ATOMS),1)
	CPPFLAGS_RAND = -DRANDOMIZE_ON_ATOMS
endif
CPPFLAGS = -Wall -I$(CBRAID_INCDIR) -L$(CBRAID_LIBDIR) $(CPPFLAGS_CLN) $(CPPFLAGS_RAND) $(CPPFLAGS_PAR) -g
LIBFLAGS = -lcbraid $(LIBFLAGS_CLN) -ltbb
MAKEDEPFLAGS = -M

# Default arguments for test programs.
TEST_ARG =

# Default compile rules
# Rule to make an object file from a source file.
%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) -I$(CBRAID_INCDIR) $<

# Rule to make an executable from an object file.
.o:
	$(CPP) $(CPPFLAGS) $(OBJ) $@.o -o $(subst _main,,$@) $(LIBFLAGS)

# Default target.
all: $(PROGRAMS)

braiding: braiding_main

$(PROGRAMS): $(DEPFILE) $(CBRAID_LIB)

# Build library
$(CBRAID_LIB):
	cd $(CBRAID_LIBDIR); make libcbraid.a

# Sample programs.
do-test: test
	./test $(TEST_ARG)

do-speedtest: speedtest
	./speedtest $(SPEEDTEST_ARG)

# Cleanup.
clean:
	rm -rf $(subst _main,,$(PROGRAMS)) $(OBJS) $(DEPFILE)
	cd $(CBRAID_LIBDIR); make clean

# Make dependencies.
depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

# Make dependencies file.
$(DEPFILE):
	$(MAKEDEPCMD) $(MAKEDEPFLAGS) $(SOURCES) \
		-I$(CBRAID_INCDIR) > $(DEPFILE)

# Read in dependencies file created by 'make depend'.
sinclude $(DEPFILE)
