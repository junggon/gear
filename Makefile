.PHONY: all clean gear debug install_dependencies_ubuntu
 
all: gear

clean: 
	cd build_release;make clean;rm CMakeCache.txt;rm Makefile

gear: build_release/Makefile
	cd build_release;make --no-print-directory

debug: build_debug/Makefile
	cd build_debug;make --no-print-directory

build_release/Makefile: CMakeLists.txt
	mkdir -p build_release;cd build_release;export CC=gcc;export CXX=g++;cmake -DCMAKE_BUILD_TYPE=Release ..

build_debug/Makefile: CMakeLists.txt
	mkdir -p build_debug;cd build_debug;export CC=gcc;export CXX=g++;cmake -DCMAKE_BUILD_TYPE=Debug ..

install_dependencies_ubuntu:
	sudo apt-get install build-essential cmake freeglut3-dev
