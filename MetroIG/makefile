CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = MetroIG
LIBS   = -lm -lgsl -lgslcblas -lpthread -lqutils -L../Lib 
OBJS   = MetroIG.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : MetroIG.h MetroIG.c
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o MetroIG
