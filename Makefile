CC=g++

CFLAGS = -c -g -std=c++0x

SOURCES = main.cpp Caller.cpp Batch.cpp Location.cpp Sample.cpp ReadcountEntry.cpp Allele.cpp Filter.cpp Statistics.cpp Common.cpp

OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLE = sinvict

all: $(SOURCES) $(EXECUTABLE)
	rm -fr *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o sinvict

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -fr $(EXECUTABLE) *.o
