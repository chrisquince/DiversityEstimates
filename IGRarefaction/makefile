CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = IGRarefaction
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = IGRarefaction.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : IGRarefaction.c IGRarefaction.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o IGRarefaction
