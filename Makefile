# This section was pulled from William Brent's convolve~ Makefile. It can be found here:
# https://github.com/wbrent/convolve_tilde/blob/main/Makefile
ldlibs = -L/usr/local/lib -lfftw3f-3 # Windows
# ldlibs = -L/usr/local/lib -lfftw3f # MacOSX / Linux
cflags = -Iinclude -I/usr/local/include # location of FFTW.h

# General purpose Makefile for Pd externals
directory = $(shell basename "$(CURDIR)")

lib.name = $(directory)
class.sources = src/$(directory).c
datafiles = $(directory)-help.pd

PDLIBBUILDER_DIR=libs/pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder