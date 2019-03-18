CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = SIAbundance
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = SIAbundance.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SIAbundance.c SIAbundance.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SIAbundance
