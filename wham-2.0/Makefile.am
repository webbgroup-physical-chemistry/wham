AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src tools 

clean: 
	-rm -rf src/*.o src/wham
	-rm -rf tools/*.o tools/bw_wham tools/pca_wham

distclean: clean
	-rm -rf bin/wham
	-rm -rf bin/bw_wham
	-rm -rf bin/pca_wham
	-rm -rf autom4te.cache
	-rm -rf tests/schemes/*/dat_*
	-rm -rf tests/schemes/*/*.h5 tests/schemes/*/*.pdf
	-rm -rf tests/schemes/*/*_list.file tests/schemes/*/*.pyc
	-rm -rf tests/include/*.pyc

install:
	mkdir -p bin
	-cp src/wham bin/wham
	-cp tools/bw_wham bin/bw_wham
	-cp tools/pca_wham bin/pca_wham

	
print-%:
	@echo '$*$($*)'
