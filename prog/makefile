CC=gcc

OBJS= \
	myprogram.o \
	libmesh5.o \
	groupfunctions.o

all: myprogram

myprogram: $(OBJS)
	$(CC) -O3 -o myprogram $(OBJS) 

libmesh5: libmesh.c
	$(CC) -c libmesh5.c -o libmesh.o

groupfunctions: groupfunctions.c libmesh.c
	$(CC) -c groupfunctions.c libmesh.c -o groupfunctions.o

clean:
	$(RM) -r $(OBJS) myprogram

