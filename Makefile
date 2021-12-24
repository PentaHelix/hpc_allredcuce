benchmark-hydra: build_benchmark
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/benchmark1 1_20x1.log
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/benchmark1 1_20x16.log
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/benchmark1 1_20x32.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=1 ./bin/benchmark1 1_32x1.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=16 ./bin/benchmark1 1_32x16.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=32 ./bin/benchmark1 1_32x32.log

	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/benchmark2 2_20x1.log
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/benchmark2 2_20x16.log
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/benchmark2 2_20x32.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=1 ./bin/benchmark2 2_32x1.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=16 ./bin/benchmark2 2_32x16.log
	srun -p q_student --time=1:00 -N 32 --ntasks-per-node=32 ./bin/benchmark2 2_32x32.log

benchmark: build_benchmark
	mpirun ./bin/benchmark1 test1.log
	mpirun ./bin/benchmark2 test2.log

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