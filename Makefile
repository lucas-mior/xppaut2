MAJOR_VERSION = 8.0
MINOR_VERSION = 1

PREFIX ?= /usr/local

CC = clang

CFLAGS = -Wall -Wextra -Wpedantic -Wfatal-errors
# CFLAGS += -Werror
# CFLAGS += -Wno-unused-parameter -Wno-unused-variable
# CFLAGS += -Wno-unused-but-set-variable -Wno-unused-but-set-parameter
# CFLAGS += -Wno-tautological-compare
# CFLAGS += -Wno-sign-compare
CFLAGS += -g -O2
CFLAGS += -I./src/ -I./bitmaps/ -I/opt/X11/include -I./
CFLAGS += -DMYSTR1=$(MAJOR_VERSION) -DMYSTR2=$(MINOR_VERSION)
CFLAGS += -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES -DHAVEDLL
CFLAGS += -D_DEFAULT_SOURCE -std=c99

LDFLAGS = -lX11 -lm -ldl -L.

SOURCES = $(wildcard src/*.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = xppaut

all: $(TARGET)

$(TARGET): $(OBJECTS) Makefile
	$(CC) $(CFLAGS) -o $(TARGET) $(filter-out Makefile, $^) $(LDFLAGS) 

%.o: %.c %.h Makefile
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.c Makefile
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o src/*.o $(TARGET)

install:
	install -Dm755 xppaut   ${DESTDIR}${PREFIX}/bin/xppaut
	install -Dm644 xppaut.1 ${DESTDIR}${PREFIX}/man/man1/xppaut.1
