run1: task1
	mpirun -np 4 task1

run2: task2
	mpirun -np 4 task2

task1: task1.c
	mpicc task1.c -lm -o task1

task2: task2.c
	mpicc task2.c -lm -o task2

clean:
	rm task1