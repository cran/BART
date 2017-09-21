
CXX = g++
CXXFLAGS = -std=c++11 -fpic -O2 -g -mtune=native
LIB = -L/usr/local/lib -lRmath 

## Rmath library for random number generation and other needs
CPPFLAGS = -I. -I/usr/local/include -DMATHLIB_STANDALONE -DRNG_Rmath

## STL random class for random number generation and the Rmath functions for other needs
##CPPFLAGS = -I. -I/usr/local/include -DMATHLIB_STANDALONE -DRNG_random

.cpp.o : Makefile common.h rn.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $*.o

pmain.out : Makefile cpbart.cpp pmain.o rtnorm.o bart.o bartfuns.o bd.o tree.o treefuns.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) pmain.o bart.o bartfuns.o bd.o tree.o treefuns.o rtnorm.o -o pmain.out $(LIB)

lmain.out : Makefile clbart.cpp lmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) lmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o -o lmain.out $(LIB)

clean :
	rm *.o *main.out
