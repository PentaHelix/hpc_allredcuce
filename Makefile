benchmark-hydra: build_benchmark
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/benchmark1
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/benchmark1
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/benchmark1
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=1 ./bin/benchmark1
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=16 ./bin/benchmark1
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=32 ./bin/benchmark1

	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/benchmark2
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/benchmark2
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/benchmark2
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=1 ./bin/benchmark2
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=16 ./bin/benchmark2
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=32 ./bin/benchmark2

benchmark: build_benchmark
	mpirun ./bin/benchmark1
	mpirun ./bin/benchmark2

validate: build_validate
	mpirun ./bin/validate1
	mpirun ./bin/validate2

build_validate: task1.o task2.o
	mpicc validate.c bin/task1.o -lm -o bin/validate1
	mpicc validate.c bin/task2.o -lm -o bin/validate2

build_benchmark: task1.o task2.o
	mpicc benchmark.c bin/task1.o -lm -o bin/benchmark1
	mpicc benchmark.c bin/task2.o -lm -o bin/benchmark2

task1.o: task1.c
	mpicc -c task1.c -lm -o bin/task1.o

task2.o: task2.c
	mpicc -c task2.c -lm -o bin/task2.o

clean:
	rm bin/*