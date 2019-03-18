CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = MetroLogStudent
LIBS   = -lm -lgsl -lgslcblas -lpthread -lqutils -L../Lib 
OBJS   = MetroLogStudent.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : MetroLogStudent.h MetroLogStudent.c
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o MetroLogStudent
