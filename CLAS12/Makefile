# Define paths relative to CLAS12 directory
SRC_DIR = analysis
BIN_DIR = executables
INCLUDE_DIR = ../include
LIB_DIR = ../lib

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Wall -I$(INCLUDE_DIR) $(shell root-config --cflags)
LDFLAGS = -L$(LIB_DIR) $(shell root-config --glibs)

# Only set these if TARGET is provided
ifneq ($(TARGET),)
    SRC_FILE = $(SRC_DIR)/$(TARGET).cpp
    OBJ_FILE = $(SRC_FILE:.cpp=.o)
    BIN_FILE = $(BIN_DIR)/$(TARGET)
    LIB_OBJ = $(LIB_DIR)/MomCorrParticle.o
endif

# Ensure a target is provided (ONLY for compilation, NOT for clean)
ifneq ($(MAKECMDGOALS), clean)
    ifndef TARGET
        $(error "Usage: make TARGET=<filename_without_extension>")
    endif
endif

# Rule to compile the MomCorrParticle class
$(LIB_DIR)/MomCorrParticle.o: $(LIB_DIR)/MomCorrParticle.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compilation rule for the main executable
$(BIN_FILE): $(SRC_FILE) $(LIB_OBJ)
	mkdir -p $(BIN_DIR)  # Ensure bin directory exists
	$(CXX) $(CXXFLAGS) $(SRC_FILE) $(LIB_OBJ) -o $(BIN_FILE) $(LDFLAGS)

# Cleanup rule (independent of TARGET)
clean:
	rm -f $(BIN_DIR)/* $(LIB_DIR)/*.o

