CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = MetroSichel
LIBS   = -lm -lgsl -lgslcblas -lpthread -lqutils -L../Lib 
OBJS   = MetroSichel.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : MetroSichel.h MetroSichel.c
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o MetroSichel
