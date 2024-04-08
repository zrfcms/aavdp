inc=-I ./include/lapacke:./include/hdf5:./include/libpng:./EBSD
lib=./lib/liblapacke.a ./lib/liblapack.a \
./lib/libcblas.a ./lib/librefblas.a ./lib/libm.a ./lib/libgfortran.dll.a ./lib/libhdf5.dll.a \
./lib/libpng.a ./lib/libz.a
CXX=g++.exe
CC=gcc.exe
src=$(wildcard ./*/*.cpp *.cpp)
obj=$(patsubst %.cpp, %.o, $(src))

main: $(obj)
	@echo "Start..."
	$(CXX) $(obj) -o $@ \
	-g $(inc) $(lib)
	@echo "End."

%.o: %.cpp
	@echo $<
	$(CXX) -MMD $< -c -o $@ -g

.PHONY: clean
clean: 
	rm $(obj)