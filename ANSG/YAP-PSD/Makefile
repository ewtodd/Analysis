# ROOT Analysis Framework Makefile

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -fPIC
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
LIBDIR = lib

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)

# Library name
LIBNAME = libWaveformAnalysis
SHAREDLIB = $(LIBDIR)/$(LIBNAME).so

# Default target
all: $(SHAREDLIB)

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# Compile object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -I$(INCDIR) -c $< -o $@

# Create shared library
$(SHAREDLIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS)

# CRITICAL: Run target depends on the library being built first
run-initial: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x InitialProcessing.cpp");'

run-calibrate: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Calibration.cpp");'

run-background: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Background.cpp");'

run-plots: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Plotting.cpp");'

run-PSD: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x PSD.cpp");'

run-optimize: $(SHAREDLIB)
	@echo "Setting up ROOT environment..."
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x OptimizeGates.cpp");'


# Clean build files
clean:
	rm -rf $(OBJDIR) $(LIBDIR)

.PHONY: all clean run run-alt

