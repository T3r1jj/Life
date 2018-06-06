life:
		mpiCC -O3 -fno-omit-frame-pointer -fopenmp life.cpp -o life

clean:
		rm -f life.o life
		
