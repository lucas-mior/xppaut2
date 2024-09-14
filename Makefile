MAJOR_VERSION = 8.0
MINOR_VERSION = 1

PREFIX ?= /usr/local

CC = clang

CFLAGS = -Wall -Wextra -Wpedantic -Wfatal-errors
CFLAGS += -Wno-unused-parameter
# CFLAGS += -Wno-unused-but-set-variable
CFLAGS += -Wno-unused-but-set-parameter
# CFLAGS += -g -O2
CFLAGS += -I./src/ -I./bitmaps/ -I/opt/X11/include -I./ -I./src/cvode/
CFLAGS += -DMYSTR1=$(MAJOR_VERSION) -DMYSTR2=$(MINOR_VERSION)
CFLAGS += -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES -DHAVEDLL
CFLAGS += -D_DEFAULT_SOURCE -std=c99

LDFLAGS = -lX11 -lm -ldl -L.

SOURCES = $(wildcard src/*.c src/cvode/*.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = xppaut

all: $(TARGET)

test: CFLAGS += -Wno-error
test: all

CFLAGS += -Werror

$(TARGET): $(OBJECTS) Makefile
	$(CC) $(CFLAGS) -o $(TARGET) $(filter-out Makefile, $^) $(LDFLAGS) 

%.o: %.c %.h Makefile
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.c Makefile
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o src/*.o src/cvode/*.o src/sbml/*.o $(TARGET)

install:
	install -Dm755 xppaut   ${DESTDIR}${PREFIX}/bin/xppaut
	install -Dm644 xppaut.1 ${DESTDIR}${PREFIX}/man/man1/xppaut.1
