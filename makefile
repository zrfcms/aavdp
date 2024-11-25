CXX=g++
SYS=win
INC=-I ./src/include/$(SYS):./src/QSPG/spglib:./src/QB:./src/QSPG:./src/MATH:./src/MODEL:./src/XRD:./src/KED:./src/DED:./src/KKD:./src/RDF:./src/SSF
LIB=./src/lib/$(SYS)/liblapacke.a ./src/lib/$(SYS)/liblapack.a ./src/lib/$(SYS)/libgfortran.dll.a ./src/lib/$(SYS)/libcblas.a ./src/lib/$(SYS)/librefblas.a ./src/lib/$(SYS)/libpng16.a ./src/lib/$(SYS)/libz.a

SRC=$(wildcard ./src/*/*/*.cpp ./src/*/*.cpp ./src/*.cpp)
OBJ=$(patsubst %.cpp, %.o, $(SRC))

./AAVDP: $(OBJ)
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