all: call_fgb

CFLAGS = -g
#CFLAGS = -g
CC = gcc $(CFLAGS)
CPP = g++ $(CFLAGS)

LIBDIR = ../call_FGb/nv/maple/C/x64
INCDIR1 = ../call_FGb/nv/maple/C
INCDIR2 = ../call_FGb/nv/protocol
INCDIR3 = ../call_FGb/nv/int

call_fgb.o: call_fgb.c
	$(CC) -I $(INCDIR1) -I $(INCDIR2) -I $(INCDIR3)  -c $<

call_fgb: call_fgb.o
	$(CPP)  -o $@ call_fgb.o  -L$(LIBDIR) -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp  -lm -fopenmp

