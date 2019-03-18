CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = LSRarefaction
LIBS   = -lm -lgsl -lgslcblas -lpthread
OBJS   = LSRarefaction.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : LSRarefaction.h LSRarefaction.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -rf *.o LSRarefaction

