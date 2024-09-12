MAJOR_VERSION = 8.0
MINOR_VERSION = 1

CC = gcc
CFLAGS = -Wall -Wextra
CFLAGS += -DMYSTR1=$(MAJOR_VERSION) -DMYSTR2=$(MINOR_VERSION)
CFLAGS += -g -pedantic -O2
CFLAGS += -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES -DHAVEDLL
CFLAGS += -D_DEFAULT_SOURCE -std=c99
CFLAGS += -I./src/ -I./bitmaps/ -I/opt/X11/include
LDFLAGS = -lX11 -lm -ldl
SOURCES = $(wildcard src/*.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = xppaut

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f src/*.o $(TARGET)
