CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = Subsample
LIBS   = -lm -lgsl -lgslcblas
OBJS   = Subsample.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : Subsample.c Subsample.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o Subsample
