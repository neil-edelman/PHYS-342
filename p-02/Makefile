PROJ := ultimatedoom
PCPR := $(PROJ).exe
OBJS := Open.o Simulation.o
SRCS := $(OBJS:.o=.c)
H    := $(OBJS:.o=.h)

CC   := gcc
CF   := -Wall -O3 -fasm -fomit-frame-pointer -ffast-math -funroll-loops -fasm -fomit-frame-pointer -ffast-math -funroll-loops -pedantic -framework OpenGL -framework GLUT
OF   := -Wall -O3 -fasm -fomit-frame-pointer -ffast-math -funroll-loops -fasm -fomit-frame-pointer -ffast-math -funroll-loops -pedantic
PC   := /usr/local/i386-mingw32-4.3.0/bin/i386-mingw32-gcc
PF   := -Wall -O3 -fasm -fomit-frame-pointer -ffast-math -funroll-loops -fasm -fomit-frame-pointer -ffast-math -funroll-loops -pedantic -DGL -lglut32 -lglu32 -lopengl32 -mwindows
LIBPC:= /usr/local/i386-mingw32-4.3.0/i386-mingw32/lib/libglut32.a /usr/local/i386-mingw32-4.3.0/i386-mingw32/lib/libglu32.a /usr/local/i386-mingw32-4.3.0/i386-mingw32/lib/libopengl32.a

default: $(PROJ)

$(PROJ): $(OBJS)
	$(CC) $(CF) -o $@ $(OBJS)

.c.o:
	$(CC) $(OF) -c $<

.PHONY: clean
clean:
	-rm $(OBJS)

.PHONY: pc
$(PC): $(OBJS)
	$(PC) $(PF) -o $@ $(SRCS) $(LIBPC)
