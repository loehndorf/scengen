CC=gcc 
CFLAGS=-Wall -DNO_DLL_DEFS
LDLIBS=-lm

OBJECTS=HKW_cubic.o HKW_sg.o matrix.o moments.o sg_HKW.o sg_functions.o

all: scengen_HKW

scengen_HKW: $(OBJECTS)
	$(LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f scengen_HKW $(OBJECTS)

