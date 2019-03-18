CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = ERarefaction
LIBS   = -lm -lgsl -lgslcblas
OBJS   = ERarefaction.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : ERarefaction.c ERarefaction.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o ERarefaction

