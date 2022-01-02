benchmark-hydra: build-benchmark
	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=1 ./bin/benchmark1 logs/1_20x1.log
	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=16 ./bin/benchmark1 logs/1_20x16.log
	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=32 ./bin/benchmark1 logs/1_20x32.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=1 ./bin/benchmark1 logs/1_36x1.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=16 ./bin/benchmark1 logs/1_36x16.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=32 ./bin/benchmark1 logs/1_36x32.log

	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=1 ./bin/benchmark2 logs/2_20x1.log
	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=16 ./bin/benchmark2 logs/2_20x16.log
	srun -p q_student --time=5:00 -N 20 --ntasks-per-node=32 ./bin/benchmark2 logs/2_20x32.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=1 ./bin/benchmark2 logs/2_36x1.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=16 ./bin/benchmark2 logs/2_36x16.log
	srun -p q_student --time=5:00 -N 36 --ntasks-per-node=32 ./bin/benchmark2 logs/2_36x32.log

validate-hydra: build-validate
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/validate1
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/validate1
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/validate1
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=1 ./bin/validate1
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=16 ./bin/validate1
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=32 ./bin/validate1

	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=1 ./bin/validate2
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=16 ./bin/validate2
	srun -p q_student --time=1:00 -N 20 --ntasks-per-node=32 ./bin/validate2
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=1 ./bin/validate2
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=16 ./bin/validate2
	srun -p q_student --time=1:00 -N 36 --ntasks-per-node=32 ./bin/validate2

profile: build-profile
	./mpip/bin/mpirun-mpip ./bin/profile1
	./mpip/bin/mpirun-mpip ./bin/profile2

benchmark: build-benchmark
	mpirun ./bin/benchmark1 ./report/logs/test1.log
	mpirun ./bin/benchmark2 ./report/logs/test2.log

validate: build-validate
	mpirun ./bin/validate1
	mpirun ./bin/validate2

build-profile: task1.o task2.o
	mpicc profile.c bin/task1.o -lm -o bin/profile1
	mpicc profile.c bin/task2.o -lm -o bin/profile2

build-validate: task1.o task2.o
	mpicc validate.c bin/task1.o -lm -o bin/validate1
	mpicc validate.c bin/task2.o -lm -o bin/validate2

build-benchmark: task1.o task2.o
	mpicc benchmark.c bin/task1.o -lm -o bin/benchmark1 -O3
	mpicc benchmark.c bin/task2.o -lm -o bin/benchmark2 -O3

run1: task1.o
	mpicc run.c bin/task1.o -lm -o bin/run1 -O3
	mpirun ./bin/run1

run2: task2.o
	mpicc run.c bin/task2.o -lm -o bin/run2 -O3
	mpirun ./bin/run2

task1.o: task1.c
	mpicc -c task1.c -o bin/task1.o -O3

task2.o: task2.c
	mpicc -c task2.c -o bin/task2.o -O3

clean:
	rm bin/*