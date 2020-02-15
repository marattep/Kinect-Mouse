
LIB = -lglut -lGLU -lfreenect -lXtst
CFLAGS=-fPIC -g -Wall `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`
INC = -I/usr/local/include/libfreenect/

kmouse_mm.out : kinect_mouse_mm.c
	gcc $(LIB) $(CFLAGS) $(INC) kinect_mouse_mm.c -o kmouse_mm.out $(LIBS)
	
clean :
	rm *.out
