CXX = g++

CFLAGS := -Wextra -Wall -O3 -mtune=native -std=c++14

INC = -I$(HOME)/include -I/opt/homebrew/opt/libomp/include

LIBS = -L$(HOME)/lib -lglsc3d_3 -lmt19937ar -L/opt/local/lib -lsdl2 -lfreetype -lpng -lc++ -L/$(HOME)/lib -framework OpenGL

LD = $(CXX)

PROGRAM = run
OBJECT = $(PROGRAM).o
SOURCE = main.cpp

$(PROGRAM): $(OBJGROUP)
	$(CC) $(CFLAGS) $(INC) $(LIBS) -o $(PROGRAM) $(SOURCE)

clean:;		@rm -f $(OBJS) $(PROGRAM) core
