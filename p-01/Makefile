PROJ := Simulation
PCPR := $(PROJ).exe
SRCS := Open.c Simulation.c
CC   := gcc
PC   := /usr/local/i386-mingw32-4.3.0/bin/i386-mingw32-gcc
CF   := -Wall -O3 -fasm -fomit-frame-pointer -ffast-math -funroll-loops -fasm -fomit-frame-pointer -ffast-math -funroll-loops -pedantic -framework OpenGL -framework GLUT
PF   := -Wall -O3 -fasm -fomit-frame-pointer -ffast-math -funroll-loops -fasm -fomit-frame-pointer -ffast-math -funroll-loops -pedantic -DGL -lglut32 -lglu32 -lopengl32 -mwindows
default:
	$(CC) $(CF) -o $(PROJ) $(SRCS)

pc:
	$(PC) $(PF) -o $(PCPR) $(SRCS)
