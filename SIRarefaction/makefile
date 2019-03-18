CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = SIRarefaction
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = SIRarefaction.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SIRarefaction.c SIRarefaction.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SIRarefaction
