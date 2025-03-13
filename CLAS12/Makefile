# Define paths
SRC_DIR = analysis
BIN_DIR = executables
INCLUDE_DIR = ../include
LIB_DIR = ../lib

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Wall -I$(INCLUDE_DIR) $(shell root-config --cflags)
LDFLAGS = -L$(LIB_DIR) $(shell root-config --glibs)

# Default target
TARGET ?=

# Ensure a target is provided
ifeq ($(TARGET),)
$(error "Usage: make TARGET=<filename_without_extension>")
endif

# Source and object files
SRC_FILE = $(SRC_DIR)/$(TARGET).cpp
OBJ_FILE = $(LIB_DIR)/MomCorrParticle.o
BIN_FILE = $(BIN_DIR)/$(TARGET)

# Rule to compile MomCorrParticle.o
$(OBJ_FILE): $(LIB_DIR)/MomCorrParticle.cpp $(INCLUDE_DIR)/MomCorrParticle.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile analysis executable
$(BIN_FILE): $(SRC_FILE) $(OBJ_FILE)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Cleanup rule
clean:
	rm -f $(BIN_DIR)/* $(LIB_DIR)/*.o

