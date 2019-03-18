CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = LSAbundance
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = LSAbundance.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : LSAbundance.h LSAbundance.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -rf *.o LSAbundance 

