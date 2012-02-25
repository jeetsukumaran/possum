AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

BUILT_SOURCES = possum_info.h

bindir = ${prefix}/bin

bin_PROGRAMS = possum@PROGRAMLABEL@

noinst_HEADERS = \
	character.hpp \
	cli.hpp \
	cmdopt.hpp \
	possum_defs.hpp \
	textutil.hpp

possum@PROGRAMLABEL@_CPPFLAGS = -I @NCL_INC_DIR@
possum@PROGRAMLABEL@_LDFLAGS = -L@NCL_LIB_DIR@/ncl
possum@PROGRAMLABEL@_LDADD = -lncl
possum@PROGRAMLABEL@_SOURCES = \
    possum.cpp \
    character.cpp \
    cmdopt.cpp \
    textutil.cpp \
	cli.cpp

possum_info.h: $(HEADERS) $(SOURCES)
	$(top_srcdir)/src/echo_possum_info.sh > possum_info.h
