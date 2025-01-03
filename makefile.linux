CXX=g++ -std=c++11
SYS=linux
INC=-I ./include:./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/MODEL:./src/XRD:./src/KED:./src/DED:./src/RDF
LIB=./lib/$(SYS)/liblapacke.a ./lib/$(SYS)/liblapack.a ./lib/$(SYS)/libcblas.a ./lib/$(SYS)/librefblas.a ./lib/$(SYS)/libgfortran.a ./lib/$(SYS)/libquadmath.a ./lib/$(SYS)/libpng16.a ./lib/$(SYS)/libz.a

SRC=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
OBJ=$(patsubst %.cpp, %.o, $(SRC))

./bin/AAVDP_linux: $(OBJ)
	@echo "Start building AAVDP on Linux platform..."
	$(CXX) $(OBJ) -o $@ \
	-g -static $(INC) $(LIB)
	@echo "End building AAVDP on Linux platform."

%.o: %.cpp
	@echo $<
	$(CXX) -MMD $< -c -o $@ -g 

.PHONY: clean
clean: 
	rm $(OBJ)
