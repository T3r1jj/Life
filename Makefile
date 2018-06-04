life:
		mpiCC -O1 -fno-omit-frame-pointer -g -fopenmp life.cpp -o life

clean:
		rm -f life.o life
		
