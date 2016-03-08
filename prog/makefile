CC=gcc

OBJS= \
	myprogram.o \
	libmesh5.o \
	groupfunctions.o \
	distanceMeshFunctions

all: myprogram

myprogram: $(OBJS)
	$(CC) -o myprogram $(OBJS) -lm 

libmesh5: libmesh.c
	$(CC) -c libmesh5.c -o libmesh.o

groupfunctions: groupfunctions.c libmesh.c
	$(CC) -c groupfunctions.c libmesh.c -o groupfunctions.o

distanceMeshFunctions: distanceMeshFunctions.c groupfunctions.c
	$(CC) -c distanceMeshFunctions.c groupfunctions.c -o distanceMeshFunctions.o

clean:
	$(RM) -r $(OBJS) myprogram

