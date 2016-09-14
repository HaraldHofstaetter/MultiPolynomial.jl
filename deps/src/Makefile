all: call_fgb # call_giac

#CFLAGS = 
CFLAGS = -g
CC = gcc $(CFLAGS)
CPP = g++ $(CFLAGS)

FGB_LIBDIR = -L../call_FGb/nv/maple/C/x64
FGB_INCDIR = -I../call_FGb/nv/maple/C \
              -I../call_FGb/nv/protocol \
              -I../call_FGb/nv/int

GIAC_LIBDIR =
GIAC_INCDIR = -I /usr/include/giac


call_fgb.o: call_fgb.c
	$(CC)  $(FGB_INCDIR)  -c $<

call_fgb: call_fgb.o
	$(CPP)  -o $@ call_fgb.o  $(FGB_LIBDIR) -lfgb -lfgbexp -lgb -lgbexp -lminpoly -lminpolyvgf -lgmp  -lm -fopenmp


call_giac.o: call_giac.C
	$(CPP)  $(GIAC_INCDIR)  -c $<

call_giac: call_giac.o
	$(CPP)  -o $@ call_giac.o  $(GIAC_LIBDIR) -lgiac -lgmp
