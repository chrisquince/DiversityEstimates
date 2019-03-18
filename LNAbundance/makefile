CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = LNAbundance
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = LNAbundance.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : LNAbundance.h LNAbundance.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -rf *.o LNAbundance 

