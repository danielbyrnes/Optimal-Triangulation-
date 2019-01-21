CC=gcc
CXX=g++ 
CPPFLAGS=-g -Wall --std=c++0x

TARGET=tool
SRCS=main.cc geometry.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: $(OBJS)
	g++ $(CPPFLAGS) -o $(TARGET) $(OBJS) 

main.o: main.cc geometry.h 
	g++ $(CPPFLAGS) -c main.cc 
  
geometry.o: geometry.cc geometry.h
	g++ $(CPPFLAGS) -c geometry.cc

clean:
	$(RM) $(TARGET) *.o 
