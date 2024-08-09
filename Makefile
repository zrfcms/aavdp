inc=-I ./src/include:./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/MODEL:./src/XRD:./src/SND:./src/SED:./src/KKD:./src/DKD:./src/RDF:./src/SSF
lib=./src/lib/liblapacke.a ./src/lib/liblapack.a ./src/lib/libcblas.a ./src/lib/librefblas.a \
./src/lib/libgfortran.dll.a ./src/lib/libpng16.a ./src/lib/libz.a
CXX=g++
CC=gcc
src=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
obj=$(patsubst %.cpp, %.o, $(src))

./AAVDP: $(obj)
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