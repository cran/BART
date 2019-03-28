
##CXX = g++
##CXXFLAGS = -std=gnu++11 -Wall -g -O2 -fpic -mtune=native
CXX = `R CMD config CXX11`
CXXFLAGS = `R CMD config CXX11STD` `R CMD config CXX11FLAGS` `R CMD config CXX11PICFLAGS` 

## Rmath library for random number generation and other needs
CPPFLAGS = `R CMD config CPPFLAGS` -I. -DMATHLIB_STANDALONE -DRNG_Rmath
LIB = `R CMD config LDFLAGS` -lRmath

## STL random class for random number generation and the Rmath functions for other needs
## CPPFLAGS = -I. -I/usr/local/include -DMATHLIB_STANDALONE -DRNG_random

.cpp.o : Makefile common.h rn.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $*.o

main.out : Makefile main.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o -o main.out $(LIB)

wmain.out : Makefile cwbart.cpp wmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) wmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o -o wmain.out $(LIB)

pmain.out : Makefile cpbart.cpp pmain.o rtnorm.o bart.o bartfuns.o bd.o tree.o treefuns.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) pmain.o bart.o bartfuns.o bd.o tree.o treefuns.o rtnorm.o -o pmain.out $(LIB)

lmain.out : Makefile clbart.cpp lmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) lmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o -o lmain.out $(LIB)

mmain.out : Makefile cmbart.cpp mmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) mmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o rand_draws.o randomkit.o latent.o -o mmain.out $(LIB)

dpmain.out : Makefile cdpmbart.cpp dpmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o dp.o dpm.o dps.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) dpmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o dp.o dpm.o dps.o -o dpmain.out $(LIB)

dpmwmain.out : Makefile cdpmwbart.cpp dpmwmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o dp.o dpm.o dps.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) dpmwmain.o bart.o bartfuns.o bd.o heterbart.o heterbartfuns.o heterbd.o tree.o treefuns.o dp.o dpm.o dps.o -o dpmwmain.out $(LIB)

clean :
	rm -f *.o *main.out

check :
	@echo CXX=$(CXX)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LIB=$(LIB)
