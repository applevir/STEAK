CC=g++ 
CFLAGS=-std=c++11 

all: steak

steak: main.o samdata.o data.o smithwaterman.o ssw.o ssw_cpp.o 
	${CC} main.o samdata.o data.o smithwaterman.o ssw.o ssw_cpp.o -o steak -fopenmp  

main.o: main.cpp ltr.hpp
	${CC} -c main.cpp -fopenmp ${CFLAGS}

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
	rm main.o samdata.o data.o smithwaterman.o ssw.o ssw_cpp.o steak
