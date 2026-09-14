# Makefile for generating R packages.
# 2011 Andrew Redd
# 2018 Andreas Karlsson
#
# Roxygen uses the roxygen2 package, and will run automatically on check and all.
#
# make clean build check

PKG_VERSION=$(shell grep -i ^version ./DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package ./DESCRIPTION | cut -d : -d \  -f 2)
R_HOME?=$(shell R RHOME)

R_FILES := $(wildcard ./R/*.R)
SRC_FILES := $(wildcard ./src/*) $(addprefix ./src/, $(COPY_SRC))
PKG_FILES := ./DESCRIPTION ./NAMESPACE $(R_FILES) $(SRC_FILES)

CPP_TEST_DIR := ./test/cpp
CPP_TEST_BUILD_DIR := $(CPP_TEST_DIR)/build
CPP_TEST_SRC := $(CPP_TEST_DIR)/callfhcrc_loop_test.cpp
CPP_TEST_BIN := $(CPP_TEST_DIR)/callfhcrc_loop_test
CPP_TEST_LOCAL_CPP_SRCS := ./src/ssim_patched.cc
CPP_TEST_LOCAL_OBJS := \
	$(CPP_TEST_BUILD_DIR)/microsimulation_patched.o \
	$(CPP_TEST_BUILD_DIR)/ssim_patched.o
CPP_TEST_INCLUDED_SRCS := \
	./src/prostata.cpp \
	./src/microsimulation_patched.h
CPP_TEST_MICROSIM_INCLUDE ?= $(shell $(R_HOME)/bin/Rscript -e 'p <- system.file("include", package = "microsimulation"); if (nzchar(p)) cat(p)')
CPP_TEST_INCLUDE := $(if $(CPP_TEST_MICROSIM_INCLUDE),-I$(CPP_TEST_MICROSIM_INCLUDE),)
CPP_TEST_DEBUGFLAGS := -g3 -O0 -fno-omit-frame-pointer -fno-inline
CPP_TEST_CXXFLAGS := -std=gnu++17 $(CPP_TEST_INCLUDE) \
	$(CPP_TEST_DEBUGFLAGS) \
	$(shell $(R_HOME)/bin/R CMD config --cppflags) \
	$(shell $(R_HOME)/bin/Rscript -e 'Rcpp:::CxxFlags()') \
	$(shell $(R_HOME)/bin/Rscript -e 'RcppArmadillo:::CxxFlags()')
CPP_TEST_LDFLAGS := \
	$(shell $(R_HOME)/bin/R CMD config --ldflags) \
	$(shell $(R_HOME)/bin/Rscript -e 'Rcpp:::LdFlags()') \
	$(shell $(R_HOME)/bin/Rscript -e 'cat(microsimulation:::LdFlags())') \
	-lgtest -lgtest_main -lpthread

.PHONY: tarball install check clean build cpp-test

tarball: $(PKG_NAME)_$(PKG_VERSION).tar.gz
$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	R CMD build .

check: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check --as-cran $(PKG_NAME)_$(PKG_VERSION).tar.gz

build: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL --build $(PKG_NAME)_$(PKG_VERSION).tar.gz

cpp-test: $(CPP_TEST_BIN)
	$(CPP_TEST_BIN)

$(CPP_TEST_BIN): $(CPP_TEST_SRC) $(CPP_TEST_INCLUDED_SRCS) $(CPP_TEST_LOCAL_OBJS)
	$(CXX) $(CPP_TEST_CXXFLAGS) -o $@ $(CPP_TEST_SRC) $(CPP_TEST_LOCAL_OBJS) $(CPP_TEST_LDFLAGS)

$(CPP_TEST_BUILD_DIR):
	mkdir -p $@

$(CPP_TEST_BUILD_DIR)/microsimulation_patched.o: ./src/microsimulation_patched.cc ./src/microsimulation_patched.h | $(CPP_TEST_BUILD_DIR)
	$(CXX) $(CPP_TEST_CXXFLAGS) -c -o $@ $<

$(CPP_TEST_BUILD_DIR)/ssim_patched.o: ./src/ssim_patched.cc ./src/ssim_patched.h | $(CPP_TEST_BUILD_DIR)
	$(CXX) $(CPP_TEST_CXXFLAGS) -c -o $@ $<

install: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

./NAMESPACE: $(R_FILES)
	Rscript -e "library(roxygen2);roxygenize()"

clean:
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -r -f $(PKG_NAME).Rcheck
	-rm -r -f ./man/*
	-rm -r -f ./NAMESPACE
	-rm -f src/*.o src/*.so
	-rm -f $(CPP_TEST_BIN)
	-rm -r -f $(CPP_TEST_BUILD_DIR)

.NOTPARALLEL: # Force disabling of -j flag
all: clean build check

.PHONY: list
list:
	@echo "R files:"
	@echo $(R_FILES)
	@echo "Source files:"
	@echo $(SRC_FILES)
