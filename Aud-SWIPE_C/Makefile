target=aswipep
prefix=/usr/local

all: aswipep

aswipep: swipe.c AudSys.c vector.c fir2.c ERBFilters.c resample.c 
	$(CC) $(CFLAGS) -o $(target) swipe.c AudSys.c vector.c fir2.c ERBFilters.c resample.c  -lm -lc -lblas -llapack -lfftw3 -lsndfile -fopenmp 	

install: aswipep
	install aswipep $(prefix)/bin

clean: aswipep
	rm $(target)
