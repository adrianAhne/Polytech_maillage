CC=gcc
CFLAG = -lm
OBJS= \
	myprogram.o \
	libmesh5.o \
	courbure.o \
	translation.o \
	superposition.o \
	rotation.o \
	normale.o \
	distanceMeshFunctions.o \
	bucket.o \
	hash.o


all: myprogram

myprogram: $(OBJS)
	$(CC) -o main $(OBJS) $(CFLAG)

libmesh5: libmesh.c
	$(CC) -c libmesh5.c -o libmesh.o


courbure: courbure.c libmesh.c
	$(CC) -c courbure.c libmesh.c -o courbure.o $(CFLAG)

translation: translation.c libmesh.c
	$(CC) -c translation.c libmesh.c -o translation.o $(CFLAG)

superposition: superposition.c libmesh.c
	$(CC) -c superposition.c libmesh.c -o superposition.o $(CFLAG)

rotation: rotation.c libmesh.c
	$(CC) -c rotation.c libmesh.c -o rotation.o $(CFLAG)

normale: normale.c libmesh.c
	$(CC) -c normale.c libmesh.c -o normale.o $(CFLAG)

distanceMeshFunctions: distanceMeshFunctions.c libmesh.c
	$(CC) -c distanceMeshFunctions.c libmesh.c -o distanceMeshFunctions.o $(CFLAG)
		


bucket: bucket.c 
	$(CC) -c bucket.c -o bucket.o $(CFLAG) 

hash: hash.c 
	$(CC) -c hash.c -o bucket.o -lm

clean:
	$(RM) -r $(OBJS) myprogram

