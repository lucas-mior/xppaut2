MAJOR_VERSION = 8.0
MINOR_VERSION = 1

CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic
CFLAGS += -DMYSTR1=$(MAJOR_VERSION) -DMYSTR2=$(MINOR_VERSION)
CFLAGS += -g -O2
CFLAGS += -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES -DHAVEDLL
CFLAGS += -D_DEFAULT_SOURCE -std=c99
CFLAGS += -I./src/ -I./bitmaps/ -I/opt/X11/include -I./
LDFLAGS = -lX11 -lm -ldl -L.
SOURCES = $(wildcard src/*.c)
OBJECTS = $(SOURCES:.c=.o)
TARGET = xppaut

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS) 

%.o: %.c %.h
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o src/*.o $(TARGET)
