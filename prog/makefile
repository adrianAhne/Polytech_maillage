CC=gcc

OBJS= \
	myprogram.o \
	libmesh5.o \
	groupfunctions.o \
	bucket.o

all: myprogram

myprogram: $(OBJS)
	$(CC) -o myprogram $(OBJS) -lm 

libmesh5: libmesh.c
	$(CC) -c libmesh5.c -o libmesh.o

groupfunctions: groupfunctions.c libmesh.c
	$(CC) -c groupfunctions.c libmesh.c -o groupfunctions.o
	
bucket: bucket.c 
	$(CC) -c bucket.c -o bucket.o -lm 
clean:
	$(RM) -r $(OBJS) myprogram

