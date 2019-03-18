CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = LNRarefaction
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = LNRarefaction.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : LNRarefaction.h LNRarefaction.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -rf *.o LNRarefaction

