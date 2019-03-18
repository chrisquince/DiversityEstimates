CC     = gcc
CFLAGS = -O3
LIBS   = -lm 
OBJS   = MatrixUtils.o FileUtils.o

libqutils.a : $(OBJS)
	$(AR) rvu  libqutils.a $(OBJS)
	ranlib libqutils.a

$(OBJS) : FileUtils.h MatrixUtils.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o *.a
