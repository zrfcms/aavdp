inc=-I ./src/include:./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/MODEL:./src/XRD:./src/SND:./src/SED:./src/KKD:./src/DKD:./src/RDF:./src/SSF
lib=./src/lib/liblapacke.a ./src/lib/liblapack.a ./src/lib/libcblas.a ./src/lib/librefblas.a \
./src/lib/libgfortran.dll.a ./src/lib/libpng16.a ./src/lib/libz.a
CXX=g++
CC=gcc
src=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
obj=$(patsubst %.cpp, %.o, $(src))

./bin/AAVDP: $(obj)
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

# inc=-I ./include/lapacke:./include/hdf5:./include/libpng:./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/HDF5:./src/MODEL:./src/XRD:./src/SND:./src/SED:./src/KKD:./src/DKD:./src/RDF:./src/SSF
# lib=./lib/liblapacke.a ./lib/liblapack.a \
# ./lib/libcblas.a ./lib/librefblas.a ./lib/libm.a ./lib/libgfortran.dll.a ./lib/libhdf5.dll.a \
# ./lib/libpng16.a ./lib/libz.a
# CXX=g++
# CC=gcc
# src=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
# obj=$(patsubst %.cpp, %.o, $(src))

# ./bin/AAVDP: $(obj)
# 	@echo "Start..."
# 	$(CXX) $(obj) -o $@ \
# 	-g $(inc) $(lib)
# 	@echo "End."

# %.o: %.cpp
# 	@echo $<
# 	$(CXX) -MMD $< -c -o $@ -g 

# .PHONY: clean
# clean: 
# 	rm $(obj)