CC = gcc

OBJS= myprogram.o libmesh5.o 

all: myprogram

myprogram: $(OBJS)
	$(CC) -o myprogram $(OBJS) 
libmesh5.o: libmesh5.c libmesh5.h
	gcc -c libmesh5.c
clean:
	$(RM) -r $(OBJS) myprogram

