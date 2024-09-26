MAJOR_VERSION = 8.0
MINOR_VERSION = 3

PREFIX ?= /usr/local

CFLAGS = -D_DEFAULT_SOURCE -std=c99 -DMALLOC_DEBUG
# CFLAGS += -fsanitize=undefined
CFLAGS += -Wall -Wextra -Wpedantic -Wfatal-errors
CFLAGS += -Wfatal-errors
CFLAGS += -I./src/ -I./bitmaps/ -I./ -I./src/cvode/
CFLAGS += -DMAJOR_VERSION=$(MAJOR_VERSION) -DMINOR_VERSION=$(MINOR_VERSION)
CFLAGS += -DNON_UNIX_STDIO -DAUTO -DCVODE_YES -DHAVEDLL

LDFLAGS = -lX11 -lm -ldl -L.

SOURCES = $(wildcard src/*.c src/cvode/*.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = xppaut

all: $(TARGET)

clang: C = clang
clang: CFLAGS += -Weverything -Wno-unsafe-buffer-usage -Wno-padded
clang: CFLAGS += -Wno-float-equal
clang: CFLAGS += -Wno-format-nonliteral
clang: all

gcc: C = gcc
gcc: all

tcc: C = tcc
tcc: CFLAGS += -Wno-error
tcc: all

bear: CFLAGS += -Wno-error
bear: compile_commands.json

compile_commands.json: Makefile
	rm -f *.o src/*.o src/cvode/*.o src/sbml/*.o $(TARGET)
	bear -- make -j1 clang > compile_commands.json

tags: CFLAGS += -Wno-error
tags: $(OBJECTS)
	ctags -o tags --kinds-C=+l src/*.c src/*.h
	vtags.sed tags > .tags.vim

meta: C = clang
meta: bear tags

test: C = $(CC)
test: CFLAGS += -Wno-error
test: all
	./xppaut & sleep 0.7
	killall xppaut
	make clean

CFLAGS += -Werror

# this files contain macros that generate false warnings
src/autlib1.o: CFLAGS += -Wno-unused-but-set-variable
src/autlib2.o: CFLAGS += -Wno-unused-but-set-variable
src/autlib3.o: CFLAGS += -Wno-unused-but-set-variable
src/autlib4.o: CFLAGS += -Wno-unused-but-set-variable
src/autlib5.o: CFLAGS += -Wno-unused-but-set-variable
src/extra.o: CFLAGS += -Wno-pedantic
src/form_ode.o: CFLAGS += -Wno-pedantic
src/parserslow.o: CFLAGS += -Wno-pedantic

$(TARGET): $(OBJECTS) Makefile
	$(C) $(CFLAGS) -o $(TARGET) $(filter-out Makefile, $^) $(LDFLAGS) 

%.o: %.c %.h Makefile
	$(C) $(CFLAGS) -o $@ -c $<

%.o: %.c Makefile
	$(C) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o src/*.o src/cvode/*.o src/sbml/*.o $(TARGET)
	rm compile_commands.json tags

install:
	install -Dm755 xppaut   ${DESTDIR}${PREFIX}/bin/xppaut
	install -Dm644 xppaut.1 ${DESTDIR}${PREFIX}/man/man1/xppaut.1
