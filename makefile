OBJS = reading_input.o main.o
DEBUG = -g3
OPTFLAG = #-O3
CC = icc $(OPTFLAG)
CFLAGS = -c $(DEBUG)
LFLAGS = $(DEBUG)
MKL_LIB = /opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a /opt/intel/mkl/lib/libmkl_sequential.a
MKL_LIB += -ldl -lpthread -lm
MKL_include = -I/opt/intel/mkl/include
OPENMP =/opt/intel/compilers_and_libraries_2016.3.170/mac/compiler/lib/
LIBS_1 =  -L$(OPENMP) 
 

all :  $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o hfexe $(MKL_include) $(MKL_LIB) $(LIBS_1)

reading_input.o : tensor_type.h reading_input.cpp
	$(CC) $(CFLAGS) reading_input.cpp $(MKL_include)

main.o : tensor_type.h main.cpp
	$(CC) $(CFLAGS) main.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

clean:
	 rm *.o hfexe
