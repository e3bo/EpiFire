CC=g++ 
CFLAGS=-c -O3 -Wall -Wno-sign-compare --pedantic
#CFLAGS=-c -g
SOURCES= Network.cpp Utility.cpp 
INCLUDE= -I./
LDFLAGS= 


OBJECTS=$(SOURCES:.cpp=.o)
LIB=libsim.a

all: $(SOURCES) $(LIB)
	
$(LIB): $(OBJECTS) 
	$(AR) -rv $(LIB) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	rm $(OBJECTS) $(LIB)

