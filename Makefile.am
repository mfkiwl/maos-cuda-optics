ACLOCAL_AMFLAGS=-I m4 --install

SUBDIRS=
SUBDIRS+=sys math lib
if USE_CUDA
SUBDIRS+=cuda
endif
SUBDIRS+=maos skyc
if USE_MEX
SUBDIRS+=mex
endif
SUBDIRS+=tools test 
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 

#For documentation only. Enable USE_COMPLEX switch for complex function documentation.
MATH_SOURCE=$(srcdir)/math/mat.cpp $(srcdir)/math/cell.cpp $(srcdir)/math/matbin.cpp $(srcdir)/math/spbin.cpp $(srcdir)/math/sp.cpp
all-local:$(srcdir)/math/tmp.cpp
$(srcdir)/math/tmp.cpp:$(MATH_SOURCE)
	@echo "#define USE_COMPLEX" > $(srcdir)/math/tmp.cpp
	@cat $(MATH_SOURCE) >> $(srcdir)/math/tmp.cpp

EXTRA_DIST = doxygen.cfg docs version.h scripts config

BUILT_SOURCES=version #non existent file will always run
version:
	@echo static const char \*GIT_VERSION=\"`cd $(srcdir) && git show -s --format='%h %ci'`\"\; >$(srcdir)/version.h

.ONESHELL:
