CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = SIShift
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = SIShift.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SIShift.c SIShift.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SIShift
