CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = MetroLogNormal
LIBS   = -lm -lgsl -lgslcblas -lpthread -lqutils -L../Lib 
OBJS   = MetroLogNormal.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : MetroLogNormal.h MetroLogNormal.c
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o MetroLogNormal
