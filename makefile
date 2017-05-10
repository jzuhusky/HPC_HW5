
CC=gcc 
flags=-O3 -lrt -lm

METHODS=serial2DMultigrid

all:${METHODS}

serial2DMultigrid: serial2DMultigrid.c
	${CC} ${flags} $^ -o serial.run

clean:
	rm -f *.run
