#let's pass the compil in make argument, eg:clang
COMPIL=gcc

CC=$(COMPIL) -g
CFLAG = -lm 
OBJS= \
	main.o \
	InOut.o \
	libmesh5.o \
	curvature.o \
	rigidtransformation.o \
	normal.o \
	ball.o \
	bucket.o \
	distance.o \
	hash.o \
	hausdorff.o


all: main

main: $(OBJS)
	$(CC) -o main $(OBJS) $(CFLAG)

libmesh5: libmesh.c 
	$(CC) -c libmesh5.c -o libmesh.o

InOut: InOut.c
	$(CC) -c InOut.c -o InOut.o $(CFLAG)

curvature: curvature.c libmesh.c
	$(CC) -c curvature.c libmesh.c -o curvature.o $(CFLAG)

rigidtransformation: rigidtransformation.c libmesh.c
	$(CC) -c rigidtransformation.c libmesh.c -o rigidtransformation.o $(CFLAG)

normal: normal.c libmesh.c
	$(CC) -c normal.c libmesh.c -o normal.o $(CFLAG)
	
ball: ball.c hash.c
	$(CC) -c ball.c -o ball.o $(CFLAG)

distance: distance.c libmesh.c bucket.c
	$(CC) -c distance.c libmesh.c bucket.c -o distance.o $(CFLAG)

bucket: bucket.c 
	$(CC) -c bucket.c -o bucket.o $(CFLAG)  

hash: hash.c bucket.c ball.c
	$(CC) -c hash.c bucket.c ball.c -o hash.o $(CFLAG)


	
hausdorff: hausdorff.c hash.c
	$(CC) -c hausdorff.c hash.c -o hausdorff.o $(CFLAG)

clean:
	$(RM) -r $(OBJS) myprogram


