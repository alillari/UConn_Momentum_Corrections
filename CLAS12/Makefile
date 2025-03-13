# Define paths relative to CLAS12 directory
SRC_DIR = analysis
BIN_DIR = executables
INCLUDE_DIR = ../include
LIB_DIR = ../libs

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -I$(INCLUDE_DIR) $(shell root-config --cflags)
LDFLAGS = -L$(LIB_DIR) $(shell root-config --glibs)

# Default target
TARGET ?= 

# Find source file based on TARGET
SRC_FILE = $(SRC_DIR)/$(TARGET).cpp
OBJ_FILE = $(SRC_FILE:.cpp=.o)
BIN_FILE = $(BIN_DIR)/$(TARGET)

# Ensure a target is provided
ifeq ($(TARGET),)
$(error "Usage: make TARGET=<filename_without_extension>")
endif

# Compilation rule
$(BIN_FILE): $(SRC_FILE)
	mkdir -p $(BIN_DIR)  # Ensure bin directory exists
	$(CXX) $(CXXFLAGS) $(SRC_FILE) -o $(BIN_FILE) $(LDFLAGS)

# Cleanup rule
clean:
	rm -f $(BIN_DIR)/*

