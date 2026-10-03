bin_PROGRAMS = pyro_kramers

dSFMT = ./dSFMT/dSFMT-common.h \
				./dSFMT/dSFMT-params.h \
				./dSFMT/dSFMT-jump.h \
				./dSFMT/dSFMT-params521.h \
				./dSFMT/dSFMT-params1279.h \
				./dSFMT/dSFMT-params2203.h \
				./dSFMT/dSFMT-params4253.h \
				./dSFMT/dSFMT-params11213.h \
				./dSFMT/dSFMT-params19937.h \
				./dSFMT/dSFMT-params44497.h \
				./dSFMT/dSFMT-params86243.h \
				./dSFMT/dSFMT-params132049.h \
				./dSFMT/dSFMT-params216091.h \
				./dSFMT/dSFMT.h \
				./dSFMT/dSFMT.c	\
				./dSFMT/dSFMT-jump.c

libs = 	./lib/in_out.c \
				./lib/initialize.c \
				./lib/measure.c \
				./lib/update.c \
				./lib/utility.c

pyro_kramers_SOURCES = 	$(dSFMT) \
										$(libs) \
										./head/header.h \
										./src/main.c

CC = icc
AM_CFLAGS = -O3 -msse2 -fopenmp -DDSFMT_MEXP=19937 -DHAVE_SSE2 -fPIE
