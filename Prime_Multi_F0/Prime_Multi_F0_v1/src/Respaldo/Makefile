target=aswipep
targetCUDA=erbKernel
testTarget=aswipept
prefix=/usr/local
CC=	mpicc
NVCC= nvcc 

all: aswipep

aswipepNoCUDA:PrimeMultiF0.c AuditiveSystem.c Clocks.c Matrix.c SPUtilities.c UnitTestsMPI.c MPIMatlabCommunicator.c main.c
	$(CC) $(CFLAGS) -o $(target) PrimeMultiF0.c AuditiveSystem.c Clocks.c Matrix.c SPUtilities.c UnitTestsMPI.c  MPIMatlabCommunicator.c main.c -lm -lc -lblas -llapack -lfftw3 -lsndfile -fopenmp -lrt
	


aswipept: PrimeMultiF0.c AuditiveSystem.c Clocks.c Matrix.c SPUtilities.c UnitTestsMPI.c   MPIMatlabCommunicator.c TestFactory.c
	$(CC) $(CFLAGS) -o $(testTarget) PrimeMultiF0.c AuditiveSystem.c Clocks.c Matrix.c SPUtilities.c UnitTestsMPI.c TestFactory.c MPIMatlabCommunicator.c -lm -lc -lblas -llapack -lfftw3 -lsndfile -fopenmp -lrt



	


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
