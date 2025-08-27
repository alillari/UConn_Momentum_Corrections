# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall $(shell root-config --cflags) -I../../include
LDFLAGS  := $(shell root-config --libs)

# Paths
LIBDIR   := ../../lib
INCDIR   := ../../include
ANADIR   := ../analysis

# Library objects (extend if you add more .cpp files in lib/)
LIBOBJS  := $(LIBDIR)/MomCorrParticle.o

# All analysis scripts (.cpp) become executables
SOURCES  := $(wildcard $(ANADIR)/*.cpp)
EXES     := $(patsubst $(ANADIR)/%.cpp, %.exe, $(SOURCES))

# Default target: build everything
all: $(EXES)

# Rule to build each executable
%.exe: $(ANADIR)/%.cpp $(LIBOBJS)
	$(CXX) $(CXXFLAGS) $< $(LIBOBJS) -o $@ $(LDFLAGS)

# Rule to build library objects
$(LIBDIR)/%.o: $(LIBDIR)/%.cpp $(INCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Housekeeping
clean:
	rm -f $(EXES) $(LIBDIR)/*.o

