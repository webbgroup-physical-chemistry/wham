AUTOMAKE_OPTIONS = foreign
SUBDIRS = src 

clean: 
	-rm -rf src/*.o src/wham
	-rm -rf bin/wham
	-rm -rf tests/schemes/*/*.h5 tests/schemes/*/*.pdf

distclean: clean
	rm -rf autom4te.cache

install:
	mkdir -p bin
	-mv src/wham bin/wham

all: install
	
