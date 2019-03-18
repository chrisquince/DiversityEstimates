CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = LNShift
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = LNShift.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : LNShift.h LNShift.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -rf *.o LNShift 

