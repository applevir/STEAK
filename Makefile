CC=mpic++
CFLAGS=-std=c++11 -fopenmp -O3
LFLAGS=-lboost_program_options -lboost_thread -lboost_system -fopenmp -O3

all: steak
debug: CFLAGS += -D"DEBUG"
debug: all

steak: parallelenvironment.o inputfile.o outputfile.o process.o samdata.o data.o smithwaterman.o ssw.o ssw_cpp.o main.o
	${CC} parallelenvironment.o inputfile.o outputfile.o process.o samdata.o data.o smithwaterman.o ssw.o ssw_cpp.o main.o -o steak ${LFLAGS}

main.o: buffer.hpp main.cpp
	${CC} -c main.cpp ${CFLAGS}

parallelenvironment.o: parallelenvironment.hpp parallelenvironment.cpp
	${CC} -c parallelenvironment.cpp ${CFLAGS}

inputfile.o: inputfile.hpp inputfile.cpp
	${CC} -c inputfile.cpp ${CFLAGS}

outputfile.o: outputfile.hpp outputfile.cpp
	${CC} -c outputfile.cpp ${CFLAGS}

process.o: buffer.hpp process.hpp process.cpp
	${CC} -c process.cpp ${CFLAGS}

samdata.o: samdata.hpp samdata.cpp
	${CC} -c samdata.cpp ${CFLAGS}

data.o: data.hpp data.cpp
	${CC} -c data.cpp ${CFLAGS}

smithwaterman.o: smithwaterman.hpp smithwaterman.cpp
	${CC} -c smithwaterman.cpp ${CFLAGS}

ssw.o: ssw.h ssw.c
	${CC} -c ssw.c ${CFLAGS}

ssw_cpp.o: ssw_cpp.h ssw_cpp.cpp
	${CC} -c ssw_cpp.cpp ${CFLAGS}

clean:
	rm steak *.o
