## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=mpic++ -std=c++11 -g
#########


## Flags to give the compiler for "release mode"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -I/media/xuejian/WORK/spectra/spectra-0.5.0/include/ -I/media/xuejian/WORK/spectra/eigen-eigen-67e894c6cd8f/







obj=main.o 
main:$(obj)
	$(CCCOM) -o main $(obj)  $(LIBSPECTRA)
main.o:main.cpp master.h
	$(CCCOM) -c main.cpp $(LIBSPECTRA)
.PHONY:clean
clean:
	rm -f main $(obj)















