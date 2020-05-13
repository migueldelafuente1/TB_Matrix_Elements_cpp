# Dale un nombre a tu programa y cuando quieras ejecutarlo haz ./bin/program
PROG := matrix_elements
# Estas lineas encuentran tus .cpp y tus .h:
HDR :=$(wildcard include/*.h)
SRC :=$(wildcard src/*.cpp)
# Esta linea define los .o que va a haber en tu proyecto a partir de los .cpp
OBJ :=$(patsubst src/%.cpp, obj/%.o, $(SRC))
BIN :=bin

# Aqui definimos toda las opctiones de la compilacion:
CXX :=g++
STD :=-std=c++11
CXXFLAGS :=-g -Wall -Wextra
INCLUDE := -Iinclude

# Esto de PHONY solamente le dice que el target $(PROG) no se corresponde con un
# fichero. No tiene mucha importancia.
.PHONY: $(PROG)
# Cuando quieras recompilar todo haz make clean
.PHONY: clean

# El funcionamiento de Makefile es el siguiente:
# Cuando haces make este archivo se ejecuta y makefile tratará de conseguir hacer
# el objetivo $(PROG). Verá que este objetivo depende de $(BIN)/$(PROG) y se irá
# a ver si este esta actualizado. En ese momento verá que $(BIN)/$(PROG) depende
# de la larga lista de ficheros .o que has guardado en $(OBJ) asi que irá uno por
# uno comprobando que están todos actualizados. Cuando vea uno que no existe lo
# compilará a partir de su fichero fuente localizado en src/.
# Cuando todos los .o estén entonces procedera a ejecutar las ordenes que hay
# dentro de $(BIN)/$(PROG). Y despues termina ya.

# Cosas importantes!!!!!!!!!!!!
# - Makefile es un coñazo y necesitas poner un tabulador dentro de cada objetivo.
#   Es decir, no vale poner espacios, tiene que ser un tabulador.
# - Fijate que tengo que crear las carpetas bin y obj
# - Cada objetivo tiene que ser el nombre de un fichero (a no ser que lo pongas
#   en .PHONY) y makefile buscará ese fichero para ver si está actualizado. Si
#   no lo encuentra entonces pensará que necesita crearlo y si no has definido
#   como debe crearlo se va a cagar en ti.
# - Estoy usando una cosa que se llaman variables automaticas: $@, $^ y $<.
# - La wildcard de Makefile es % y significa: sustituye aqui el texto que quieras.

$(PROG): $(BIN)/$(PROG)

$(BIN)/$(PROG): $(OBJ)
	mkdir -p bin
	g++ $(STD) $(CXXFLAGS) $(INCLUDE) -o $@ $^

# Estoy usando una regla implicita. Es una regla general que dice a cualquier
# fichero con este patron buscalé la dependencia con ese otro patron y compilalo
# de esta manera.
obj/%.o: src/%.cpp $(HDR)
	mkdir -p obj
	g++ $(STD) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<

clean:
	rm -rf bin obj
