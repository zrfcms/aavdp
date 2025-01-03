CXX=g++
SYS=win
INC=-I ./include:./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/MODEL:./src/XRD:./src/KED:./src/DED:./src/RDF
LIB=./lib/$(SYS)/liblapacke.a ./lib/$(SYS)/liblapack.a ./lib/$(SYS)/libgfortran.dll.a ./lib/$(SYS)/libcblas.a ./lib/$(SYS)/librefblas.a ./lib/$(SYS)/libpng16.a ./lib/$(SYS)/libz.a

SRC=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
OBJ=$(patsubst %.cpp, %.o, $(SRC))

./bin/AAVDP_win: $(OBJ)
	@echo "Start building AAVDP..."
	$(CXX) $(OBJ) -o $@ \
	-g $(INC) $(LIB)
	@echo "End building AAVDP."

%.o: %.cpp
	@echo $<
	$(CXX) -MMD $< -c -o $@ -g 

.PHONY: clean
clean: 
	rm $(OBJ)