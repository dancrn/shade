ifeq ($(CXX),)
	cxx=g++
else
	cxx=$(CXX)
endif

ifeq ($(shell uname), Darwin)
	ldflags=-framework SDL2 -framework OpenGL
else
	#assuming linux..
	ldflags=-lSDL2 -lGL
endif
