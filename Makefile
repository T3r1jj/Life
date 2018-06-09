life:
		mpiCC -O3 -fopenmp life.cpp -o life

clean:
		rm -f life.o life
		
