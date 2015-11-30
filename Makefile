#CC=mpicxx
CC=mpixlcxx_r

CFLAGS_VIZ=-fopenmp -g -I/gpfs/bbp.cscs.ch/apps/viz/bbp/dev/hdf5/1.8.15/include -lz
CFLAGS_BG=-qsmp=omp:noopt -O0 -g -I/gpfs/bbp.cscs.ch/apps/bgq/external/hdf5/hdf5-1.8.5/install/include 

LIBS_VIZ=/gpfs/bbp.cscs.ch/apps/viz/bbp/dev/hdf5/1.8.15/lib/libhdf5.a
LIBS_BG=-L/gpfs/bbp.cscs.ch/apps/bgq/external/hdf5/hdf5-1.8.5/install/lib -lhdf5

viz: stats_viz.x

bg: stats_bg.x

stats_bg.x: main.cpp NESTNodeSynapse_bg.o stopwatch_bg.o
	$(CC) $(CFLAGS_BG) main.cpp -o stats_bg.x NESTNodeSynapse_bg.o stopwatch_bg.o $(LIBS_BG)

stats_viz.x: main.cpp NESTNodeSynapse_viz.o stopwatch_viz.o
	$(CC) $(CFLAGS_VIZ) main.cpp -o stats_viz.x NESTNodeSynapse_viz.o stopwatch_viz.o $(LIBS_VIZ)

NESTNodeSynapse_bg.o: NESTNodeSynapse.cpp
	$(CC) $(CFLAGS_BG) -c NESTNodeSynapse.cpp -o NESTNodeSynapse_bg.o
	
stopwatch_bg.o: timer/stopwatch.cpp
	$(CC) $(CFLAGS_BG) -c timer/stopwatch.cpp -o stopwatch_bg.o
	
	
NESTNodeSynapse_viz.o: NESTNodeSynapse.cpp
	$(CC) $(CFLAGS_VIZ) -c NESTNodeSynapse.cpp -o NESTNodeSynapse_viz.o
	
stopwatch_viz.o: timer/stopwatch.cpp
	$(CC) $(CFLAGS_VIZ) -c timer/stopwatch.cpp -o stopwatch_viz.o
	
clean: 
	rm -rf *.o *.x
