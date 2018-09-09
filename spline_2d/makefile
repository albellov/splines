CC=g++
CFLAGS=-I.

data_gen : tools\data_gen.cpp
	$(CC) -o data_gen tools\data_gen.cpp -std=c++11 $(CFLAGS)

spline : tools\tools.cpp src\spline.cpp src\main.cpp
	$(CC) -o spline src\main.cpp src\spline.cpp tools\tools.cpp -std=c++11 $(CFLAGS)

all : data_gen spline
	