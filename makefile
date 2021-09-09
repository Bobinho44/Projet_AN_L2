# Modèle de fichier Makefile, à adapter pour le TP
# Faire une copie de ce fichier, changer le nom pour "makefile" et l'adapter pour votre projet
# >cp makefile_model.txt makefile
# 
# Maintenant pour compiler il suffit d'écrire 
# >make

# options de compilation
CC = g++
CCFLAGS = -Wall -std=c++11
LIBS = -lm -lstdc++	# par exemple, -lm rajoute le libm standard
LIBSDIR = 

# fichiers du projet
SRC = matrix.cpp main.cpp remez.cpp util.cpp
OBJ = $(SRC:.c=.o)
EXEC = run.out


# règle initiale
all: $(EXEC)

# dépendance des .h
main.o: remez.hpp
matrix.o: matrix.hpp
remez.o: matrix.hpp util.hpp
util.o: util.hpp

# règles de compilation
%.o: %.c
	$(CC) $(CCFLAGS) -o $@ -c $<
	
# règles d'édition de liens
$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS) $(LIBSDIR)

# règles supplémentaires
clean:
	rm -f *.o
mproper:
	rm -f $(EXEC) *.o