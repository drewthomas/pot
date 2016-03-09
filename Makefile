MPICC=mpicc
CFLAGS=-Wextra -Wall -ansi -std=c99 -lm
MANXSLT=/usr/share/xml/docbook/stylesheet/docbook-xsl/manpages/docbook.xsl

POTSRCS=pot.c prng.c  # coord.c's included by pot.c, so don't list it here

pot-with-doc: pot pot-doc

pot: $(POTSRCS) coord.c Makefile
	$(MPICC) -s -Os $(POTSRCS) $(CFLAGS) -o pot
#	$(MPICC) -s -Os $(POTSRCS) $(CFLAGS) -lGL -lGLU -lglfw -o pot
#	$(MPICC) -g3 $(POTSRCS) $(CFLAGS) -o pot
#	$(MPICC) -g3 -pg $(POTSRCS) $(CFLAGS) -DFPE -o pot
#	$(MPICC) -g3 -pg $(POTSRCS) $(CFLAGS) -lGL -lGLU -lglfw -DFPE -o pot

pot-doc: doc/pot.1 doc/pot.ps

doc/pot.1 doc/pot.ps: doc/pot.xml
	xsltproc -o doc/pot.1 $(MANXSLT) doc/pot.xml
	groff -t -e -mandoc -Tps doc/pot.1 > doc/pot.ps
