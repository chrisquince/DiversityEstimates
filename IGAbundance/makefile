CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = IGAbundance
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = IGAbundance.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : IGAbundance.c IGAbundance.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o IGAbundance
