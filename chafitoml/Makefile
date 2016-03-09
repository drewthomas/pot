CFLAGS=-Wall -Wextra -ansi -std=c99

all: chafitsoml chafitoml

chafitsoml: chafitsoml.c Makefile
	gcc ${CFLAGS} -s -Os chafitsoml.c -lm -o chafitsoml
#	gcc -DDEBUG -s -Os ${CFLAGS} chafitsoml.c -lm -o chafitsoml

chafitoml: chafitoml.c Makefile
	gcc ${CFLAGS} -s -Os chafitoml.c -lm -o chafitoml
#	gcc '-DETA="eta.dat"' -s -Os ${CFLAGS} chafitoml.c -lm -o chafitoml
#	gcc -DDEBUG -s -Os ${CFLAGS} chafitoml.c -lm -o chafitoml
#	gcc -DDEBUG -g3 -ggdb -pg ${CFLAGS} chafitoml.c -lm -o chafitoml
