# Makefile for
# DPD simulations

DIR1 = .

st     =  $(DIR1)/stitcher.o  $(DIR1)/smalloc.o 
st2	   =  $(DIR1)/stitcher2.o  $(DIR1)/smalloc.o 

#CFLAGS = -O3  -Wall  -Wno-unused -I$(DIR2)/include 
#CFLAGS  = -debug -Wall -I$(DIR2)/include 
LIBS    = -lm -lc
CC      = gcc
##CFLAGS =  -Wno-unused


default:
	@echo Please specify what to make:
	@echo "  t - test of random number generator"

###Primaries:

stitcher:    $(st)
	     $(CC) $(CFLAGS) -o $@ $(st) $(LIBS)
stitcher2:    $(st2)
	     $(CC) $(CFLAGS) -o $@ $(st2) $(LIBS)		 
c:
	rm *.o stitcher stitcher2

