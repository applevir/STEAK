all:
	cd steak-paired && make
	cd steak-unpaired && make

clean:
	cd steak-paired && make clean
	cd steak-unpaired && make clean
