ifeq ($(default),1)
	CPU=1
	GPU=0
	OPENACC=0
endif

ifeq ($(CPU),1)
ifeq ($(OPENACC),1)
	CC = pgcc
	CFLAGS = -acc -Minfo -ta=multicore
	LDFLAGS = -lm -lgsl -lcblas
else
	CC = gcc
	CFLAGS = -fopenmp -std=c99 -I/usr/include -Ofast -Wno-unused-result
	LDFLAGS = -lm -lgsl -lcblas
endif
else
	CC = pgcc
	CFLAGS =  -acc -fast -Mautoinline -Minfo=accel -ta=tesla:cuda8.0,fastmath,maxregcount:255
	LDFLAGS = -lm -lgsl -lcblas
endif

SRCHARM=main.c core.c GRmath.c integrator.c metric.c radiative_transfer.c raptor_harm_model.c utilities.c  j_nu.c  rcarry.c newtonraphson.c
OBJHARM=main.o core.o GRmath.o integrator.o metric.o radiative_transfer.o raptor_harm_model.o utilities.o  j_nu.o  rcarry.o newtonraphson.o
HDRHARM=raptor_harm_model.h

harm: $(OBJHARM) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJHARM) $(LDFLAGS)

clean:
	rm *.o
	rm RAPTOR
