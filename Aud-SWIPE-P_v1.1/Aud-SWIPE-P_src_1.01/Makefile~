target=aswipep
targetCUDA=erbKernel
testTarget=aswipept
prefix=/usr/local
CC=	mpicc
NVCC= nvcc 

all: aswipep

aswipepNoCUDA: main.c swipe.c AudSys.c vector.c fir2.c ERBFilters.c resample.c MPI_Matlab_Communicator.c
	$(CC) $(CFLAGS) -o $(target) main.c swipe.c AudSys.c vector.c fir2.c  ERBFilters.c resample.c  MPI_Matlab_Communicator.c -lm -lc -lblas -llapack -lfftw3 -lsndfile -fopenmp 
	


aswipept: runTests.c swipe.c AudSys.c vector.c fir2.c ERBFilters.c resample.c MPI_Matlab_Communicator.c
	$(CC) $(CFLAGS) -o $(testTarget) runTests.c swipe.c AudSys.c vector.c fir2.c  ERBFilters.c resample.c  MPI_Matlab_Communicator.c -lm -lc -lblas -llapack -lfftw3 -lsndfile -fopenmp 



	


installTests: aswipept
	install aswipept $(prefix)/bin

install: aswipepNoCUDA
	install aswipep $(prefix)/bin
	
installCUDA: aswipepCUDA
	install aswipep $(prefix)/bin


clean: aswipep
	rm $(target)

cleanTests: aswipept
	rm $(testTarget)
