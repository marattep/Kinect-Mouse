/*
 * This file is part of the OpenKinect Project. http://www.openkinect.org
 *
 * Copyright (c) 2010 individual OpenKinect contributors. See the CONTRIB file
 * for details.
 *
 * This code is licensed to you under the terms of the Apache License, version
 * 2.0, or, at your option, the terms of the GNU General Public License,
 * version 2.0. See the APACHE20 and GPL2 files for the text of the licenses,
 * or the following URLs:
 * http://www.apache.org/licenses/LICENSE-2.0
 * http://www.gnu.org/licenses/gpl-2.0.txt
 *
 * If you redistribute this file in source form, modified or unmodified, you
 * may:
 *   1) Leave this header intact and distribute it under the same terms,
 *      accompanying it with the APACHE20 and GPL20 files, or
 *   2) Delete the Apache 2.0 clause and accompany it with the GPL2 file, or
 *   3) Delete the GPL v2 clause and accompany it with the APACHE20 file
 * In all cases you must keep the copyright notice intact and include a copy
 * of the CONTRIB file.
 *
 * Binary distributions must follow the binary distribution requirements of
 * either License.
 * 
 * This file has been additionally modified by: Tim Flaman
 */

 /* Dec 2019
 This is the mod to implement MagicMirror Mouse ad Swipe module

*/

/* June 2020
 This is the mod to implement Median Claculation and WebSocket Communications
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ncurses.h>
#include "libfreenect.h"

#include <assert.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/extensions/XTest.h>

#include <pthread.h>

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <math.h>
#include <time.h>
//#include <gsl/gsl_math.h>

// Median Calculation 
#define SCREEN (DefaultScreen(display))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
 
//Customize for your data Item type
typedef int Item;
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)
 
typedef struct Mediator_t
{
   Item* data;  //circular queue of values
   int*  pos;   //index into `heap` for each value
   int*  heap;  //max/median/min heap holding indexes into `data`.
   int   N;     //allocated size.
   int   idx;   //position in circular queue
   int   ct;    //count of items in queue
} Mediator;
 
/*--- Helper Functions ---*/
 
#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap 
 
int depth;
char *display_name;

typedef struct {
    int NearPixel_TooClose;
    int NearPixel_TooFarOrNoise;
    int freenect_log_level;
    int freenect_angle;
    int freenect_led;
    int gesture_click_area;
    int hovering_threshold;
    int near_threshold;
    int far_threshold;
    int MMM_Output_log;
    int MMM_Output_status;
    int MMM_Output_clicks;
    int MMM_Output_coords;
    int MMM_Output_swipes;
    int jsonout;
    int debug;
    int debugstop;
    int minimum_stroke_points;
    int maximum_stroke_points;
    long h_varmax;
    long v_varmax;
    int ScreenCenterX;
    int ScreenCenterY;
    int ShowScreen;
} KinectMouseConfig;

typedef struct {
    int depth;
    char *display_name;
    Display *display;
    Window main_window;
    Window root_window;
    pthread_t freenect_thread;
    volatile int die;
    int g_argc;
    char **g_argv;
    int window;
    float tmprot;
    pthread_mutex_t gl_backbuf_mutex;
    uint8_t gl_depth_front[640*480*4];
    uint8_t gl_depth_back[640*480*4];
    uint8_t gl_rgb_front[640*480*4];
    uint8_t gl_rgb_back[640*480*4];
    GLuint gl_depth_tex;
    GLuint gl_rgb_tex;
    freenect_context *f_ctx;
    freenect_device *f_dev;
    float pointerx, pointery;
    float mousex, mousey;
    int screenw, screenh;
    int stroke_x[1000], stroke_y[1000];
    int h_stroke_left2right_count, h_stroke_right2left_count;
    int v_stroke_up2down_count, v_stroke_down2up_count;
    long h_stroke_left2right_sum, h_stroke_right2left_sum;
    long v_stroke_up2down_sum, v_stroke_down2up_sum;
    long h_sum, v_sum;
    float h_mean, v_mean, h_variance, v_variance;
    long left2rightmul, right2leftmul, up2downmul, down2upmul;
    float DistCen[640][480];
    int current_hovering_cycles;
    int PointerX, PointerY;
    pthread_cond_t gl_frame_cond;
    int got_frames;
    int current_pixel;
    int StrokeEval;
    KinectMouseConfig config;
} KinectMouseState;

KinectMouseState state = {
    .gl_backbuf_mutex = PTHREAD_MUTEX_INITIALIZER,
    .gl_frame_cond = PTHREAD_COND_INITIALIZER,
    .die = 0,
    .tmprot = 1,
    .pointerx = 0,
    .pointery = 0,
    .mousex = 0,
    .mousey = 0,
    .screenw = 0,
    .screenh = 0,
    .got_frames = 0,
    .current_pixel = 0,
    .StrokeEval = 0,
    .config = {
        .NearPixel_TooClose = 0,
        .NearPixel_TooFarOrNoise = 0,
        .freenect_log_level = 0,
        .freenect_angle = -30,
        .freenect_led = 1,
        .gesture_click_area = 15,
        .hovering_threshold = 15,
        .near_threshold = 0,
        .far_threshold = 0,
        .MMM_Output_log = 1,
        .MMM_Output_status = 1,
        .MMM_Output_clicks = 1,
        .MMM_Output_coords = 1,
        .MMM_Output_swipes = 1,
        .jsonout = 1,
        .debug = 1,
        .debugstop = 0,
        .h_varmax = 100,
        .v_varmax = 100,
        .minimum_stroke_points = 0,
        .maximum_stroke_points = 0,
        .ScreenCenterX = 320,
        .ScreenCenterY = 240,
        .ShowScreen = 0
    }
};

/*--- Helper Functions ---*/
 
#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap 
 
int depth;
char *display_name;

Display *display;
Window main_window;
Window root_window;

pthread_t freenect_thread;
volatile int die = 0;

int g_argc;
char **g_argv;

int window;

float tmprot = 1;

pthread_mutex_t gl_backbuf_mutex = PTHREAD_MUTEX_INITIALIZER;

uint8_t gl_depth_front[640*480*4];
uint8_t gl_depth_back[640*480*4];

uint8_t gl_rgb_front[640*480*4];
uint8_t gl_rgb_back[640*480*4];

GLuint gl_depth_tex;
GLuint gl_rgb_tex;

freenect_context *f_ctx;
freenect_device *f_dev;

// The selected code is redundant and can be deleted as these variables are already part of the KinectMouseState struct.

//Median Calculation Functions

//returns 1 if heap[i] < heap[j]
int mmless(Mediator* m, int i, int j)
{
   return ItemLess(m->data[m->heap[i]],m->data[m->heap[j]]);
}
 
//swaps items i&j in heap, maintains indexes
int mmexchange(Mediator* m, int i, int j)
{
   int t = m->heap[i];
   m->heap[i]=m->heap[j];
   m->heap[j]=t;
   m->pos[m->heap[i]]=i;
   m->pos[m->heap[j]]=j;
   return 1;
}
 
//swaps items i&j if i<j;  returns true if swapped
int mmCmpExch(Mediator* m, int i, int j)
{
   return (mmless(m,i,j) && mmexchange(m,i,j));
}
 
//maintains minheap property for all items below i/2.
void minSortDown(Mediator* m, int i)
{
   for (; i <= minCt(m); i*=2)
   {  if (i>1 && i < minCt(m) && mmless(m, i+1, i)) { ++i; }
      if (!mmCmpExch(m,i,i/2)) { break; }
   }
}
 
//maintains maxheap property for all items below i/2. (negative indexes)
void maxSortDown(Mediator* m, int i)
{
   for (; i >= -maxCt(m); i*=2)
   {  if (i<-1 && i > -maxCt(m) && mmless(m, i, i-1)) { --i; }
      if (!mmCmpExch(m,i/2,i)) { break; }
   }
}
 
//maintains minheap property for all items above i, including median
//returns true if median changed
int minSortUp(Mediator* m, int i)
{
   while (i>0 && mmCmpExch(m,i,i/2)) i/=2;
   return (i==0);
}
 
//maintains maxheap property for all items above i, including median
//returns true if median changed
int maxSortUp(Mediator* m, int i)
{
   while (i<0 && mmCmpExch(m,i/2,i))  i/=2;
   return (i==0);
}
 
/*--- Public Interface ---*/
 
 
//creates new Mediator: to calculate `nItems` running median. 
//mallocs single block of memory, caller must free.
Mediator* MediatorNew(int nItems)
{
   int size = sizeof(Mediator)+nItems*(sizeof(Item)+sizeof(int)*2);
   Mediator* m=  malloc(size);
   m->data= (Item*)(m+1);
   m->pos = (int*) (m->data+nItems);
   m->heap = m->pos+nItems + (nItems/2); //points to middle of storage.
   m->N=nItems;
   m->ct = m->idx = 0;
   while (nItems--)  //set up initial heap fill pattern: median,max,min,max,...
   {  m->pos[nItems]= ((nItems+1)/2) * ((nItems&1)?-1:1);
      m->heap[m->pos[nItems]]=nItems;
   }
   return m;
}
 
 
//Inserts item, maintains median in O(lg nItems)
void MediatorInsert(Mediator* m, Item v)
{
   int isNew=(m->ct<m->N);
   int p = m->pos[m->idx];
   Item old = m->data[m->idx];
   m->data[m->idx]=v;
   m->idx = (m->idx+1) % m->N;
   m->ct+=isNew;
   if (p>0)         //new item is in minHeap
   {  if (!isNew && ItemLess(old,v)) { minSortDown(m,p*2);  }
      else if (minSortUp(m,p)) { maxSortDown(m,-1); }
   }
   else if (p<0)   //new item is in maxheap
   {  if (!isNew && ItemLess(v,old)) { maxSortDown(m,p*2); }
      else if (maxSortUp(m,p)) { minSortDown(m, 1); }
   }
   else            //new item is at median
   {  if (maxCt(m)) { maxSortDown(m,-1); }
      if (minCt(m)) { minSortDown(m, 1); }
   }
}
 
//returns median item (or average of 2 when item count is even)
Item MediatorMedian(Mediator* m)
{
   Item v= m->data[m->heap[0]];
   if ((m->ct&1)==0) { v= ItemMean(v,m->data[m->heap[-1]]); }
   return v;
}
 
 
/*--- Test Code ---*/
void PrintMaxHeap(Mediator* m)
{
   int i;
   if(maxCt(m))
      printf("Max: %3d",m->data[m->heap[-1]]);
   for (i=2;i<=maxCt(m);++i)
   {
      printf("|%3d ",m->data[m->heap[-i]]);
      if(++i<=maxCt(m)) printf("%3d",m->data[m->heap[-i]]);
   }
   printf("\n");
}
void PrintMinHeap(Mediator* m)
{
   int i;
   if(minCt(m))
      printf("Min: %3d",m->data[m->heap[1]]);
   for (i=2;i<=minCt(m);++i)
   {
      printf("|%3d ",m->data[m->heap[i]]);
      if(++i<=minCt(m)) printf("%3d",m->data[m->heap[i]]);
   }
   printf("\n");
}
 
void ShowTree(Mediator* m)
{
   PrintMaxHeap(m);
   printf("Mid: %3d\n",m->data[m->heap[0]]);
   PrintMinHeap(m);
   printf("\n");
}
 
int mediantest(int argc, char* argv[])
{
   int i,v;
   Mediator* m = MediatorNew(5);
 
   for (i=0;i<20;i++)
   {
      v = rand()&127;
//      v = i;
      printf("Inserting %3d \n",v);
      MediatorInsert(m,v);
      v=MediatorMedian(m);
      printf("Median = %3d.\n\n",v);
      ShowTree(m);
   }
}

//Kinect Functions

void DrawGLScene()
{
	pthread_mutex_lock(&state.gl_backbuf_mutex);

	while (state.got_frames < 2) {
		pthread_cond_wait(&state.gl_frame_cond, &state.gl_backbuf_mutex);
	}

	memcpy(state.gl_depth_front, state.gl_depth_back, sizeof(state.gl_depth_back));
	memcpy(state.gl_rgb_front, state.gl_rgb_back, sizeof(state.gl_rgb_back));
	state.got_frames = 0;
	pthread_mutex_unlock(&state.gl_backbuf_mutex);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glEnable(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, state.gl_depth_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 640, 480, 0, GL_RGB, GL_UNSIGNED_BYTE, state.gl_depth_front);

	glTranslated(1280, 0, 0);
	glScalef(-1, 1, 1);

	glBegin(GL_TRIANGLE_FAN);
	glColor4f(255.0f, 255.0f, 255.0f, 255.0f);
	glTexCoord2f(0, 0); glVertex3f(0,0,0);
	glTexCoord2f(1, 0); glVertex3f(640,0,0);
	glTexCoord2f(1, 1); glVertex3f(640,480,0);
	glTexCoord2f(0, 1); glVertex3f(0,480,0);
	glEnd();

	glBindTexture(GL_TEXTURE_2D, state.gl_rgb_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 640, 480, 0, GL_RGB, GL_UNSIGNED_BYTE, state.gl_rgb_front);

	glBegin(GL_TRIANGLE_FAN);
	glColor4f(255.0f, 255.0f, 255.0f, 255.0f);
	glTexCoord2f(0, 0); glVertex3f(640,0,0);
	glTexCoord2f(1, 0); glVertex3f(1280,0,0);
	glTexCoord2f(1, 1); glVertex3f(1280,480,0);
	glTexCoord2f(0, 1); glVertex3f(640,480,0);
	glEnd();

	glutSwapBuffers();
}

void keyPressed(unsigned char key, int x, int y)
{
	switch(key) {
		case 27:
			state.die = 1;
			pthread_join(state.freenect_thread, NULL);
			glutDestroyWindow(state.window);
			pthread_exit(NULL);
		break;
		case 'w':
			if (state.config.freenect_angle < 29) state.config.freenect_angle++;
			freenect_set_tilt_degs(state.f_dev,state.config.freenect_angle);
			printf("{ \"status\" : \"Angle: %d degrees\"}\n", state.config.freenect_angle);
		break;
		case 's':
			state.config.freenect_angle = 0;
			freenect_set_tilt_degs(state.f_dev,state.config.freenect_angle);
			printf("{ \"status\" : \"Angle: %d degrees\"}\n", state.config.freenect_angle);
		break;
		case 'x':
			if (state.config.freenect_angle > -30) state.config.freenect_angle--;
			freenect_set_tilt_degs(state.f_dev,state.config.freenect_angle);
			printf("{ \"status\" : \"Angle: %d degrees\"}\n", state.config.freenect_angle);
		break;
		case '1':
			freenect_set_led(state.f_dev,LED_GREEN);
			printf("{ \"status\" : \"LED Green\"}\n");
		break;
		case '2':
			freenect_set_led(state.f_dev,LED_RED);
			printf("{ \"status\" : \"LED Red\"}\n");
		break;
		case '3':
			freenect_set_led(state.f_dev,LED_YELLOW);
			printf("{ \"status\" : \"LED Yellow\"}\n");
		break;
		case '4':
			freenect_set_led(state.f_dev,LED_BLINK_YELLOW);
			printf("{ \"status\" : \"LED Blink Yellow\"}\n");
		break;
		case '5':
			freenect_set_led(state.f_dev,LED_BLINK_GREEN);
			printf("{ \"status\" : \"LED Blink Green\"}\n");
		break;
		case '6':
			freenect_set_led(state.f_dev,LED_BLINK_RED_YELLOW);
			printf("{ \"status\" : \"LED Blink Red Yellow\"}\n");
		break;
		case '0':
			freenect_set_led(state.f_dev,LED_OFF);
			printf("{ \"status\" : \"LED Off\"}\n");
		break;
		case 'o':
			state.tmprot+=0.1;
			printf("{ \"status\" : \"Rotation: %d degrees\"}\n", state.tmprot);
		break;
		case 'p':
			state.tmprot-=0.1;
			printf("{ \"status\" : \"Rotation: %d degrees\"}\n", state.tmprot);
		break;
	}

}

void ReSizeGLScene(int Width, int Height)
{
	glViewport(0,0,Width,Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho (0, 1280, 480, 0, -1.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW);
}

void InitGL(int Width, int Height)
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);
	glGenTextures(1, &state.gl_depth_tex);
	glBindTexture(GL_TEXTURE_2D, state.gl_depth_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glGenTextures(1, &state.gl_rgb_tex);
	glBindTexture(GL_TEXTURE_2D, state.gl_rgb_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	ReSizeGLScene(Width, Height);
}

/// @brief 
/// @param arg 
/// @return 
/**
 * @brief Thread function for handling OpenGL rendering.
 *
 * This function is intended to be run as a separate thread. It handles
 * the OpenGL rendering loop, updating the display based on the data
 * provided through the argument.
 *
 * @param arg A pointer to the data required for the rendering loop.
 *            The exact type and structure of this data should be
 *            documented elsewhere in the code.
 *
 * @return A pointer to the result of the thread execution. The exact
 *         type and meaning of this result should be documented
 *         elsewhere in the code.
 */
void *gl_threadfunc(void *arg)
{
	if(state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"OpenGL Window Opened\"} \n");
	if(state.config.debug) printf("OpenGL Window Opened\n");
	glutInit(&state.g_argc, state.g_argv);
	
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowPosition(-640, 0);
	
	if (!state.config.ShowScreen) 
		glutInitWindowSize(1, 1);
	else
		glutInitWindowSize(1280, 480);
	state.window = glutCreateWindow("Virtual Mouse for Magic Mirror");

	glutDisplayFunc(&DrawGLScene);
	glutIdleFunc(&DrawGLScene);
	glutReshapeFunc(&ReSizeGLScene);
	glutKeyboardFunc(&keyPressed);

	InitGL(1280, 480);

	glutMainLoop();

	return NULL;
}

uint16_t t_gamma[2048];

void depth_cb(freenect_device *dev, void *v_depth, uint32_t timestamp)
{
	
	// this is a callback function in the standard OpenKinect Framework returning a frame when ready
	// In this section we search for the currently pointed screen pixel
	
	// In this part of the function we analyze the frame returned and look for a "blob" of NearPixels
	// A NearPixel is a Pixel that is reported within the given range by Kinect
	// We search in every frame returned by the kinect for a group of NearPixels.
	// Hopefully this is the forearm of the subject facing the kinect
	// We assume the direction the mouse pointer is the furthest point of the blob from screen center
	// NOTE: Pixel variables are local to the function; Stroke variables are GLOBAL and persist across function call

	int i,j;
	long NearPixelXSum, NearPixelYSum; // Sum of X and Y coordinates of NearPixels. Used for Swipe evaluation
	int NearPixelX , NearPixelY; // NearPixelFound : near pixel coordinates
	int NearPixelMaxX  , NearPixelMaxY; // NearPixelFound near pixel MAX coordinates
	int NearPixelMinX  , NearPixelMinY; // NearPixelFound near pixel min coordinates
	int tx , ty;  // Pointer to current pixel analyzed in current frame
	uint16_t *depth = v_depth;
	float dpixel=0; 	// Distance of current pixel from the center of the screen
	float dpixelmax=-1;  //Distance of the furthest Pixel of Current blob of NP frame from screen center
	int mx , my;  // mouse x and y coordinates
	long NearPixelCount; // count of red pixels (near)
	int pval,pmax; // pval is current pixel depth according to kinect sensor. pmax is the maximun pval found
	pthread_mutex_lock(&state.gl_backbuf_mutex);

	if(state.config.debug) printf("___________________________BEGINOFRAME_________________________\n");
	if(state.config.debug) printf("Got a Frame, Anlyzing it\n");
//
//Loop all pixels of current frame and search for at least one near pixel
//
	NearPixelCount=0;
	tx=0;
	ty=0;
	NearPixelXSum=0;
	NearPixelYSum=0;
	NearPixelX = 0;
	NearPixelY = 0; 
	NearPixelMaxX = -2000;
	NearPixelMaxY = -2000;
	NearPixelMinX = 2000;
	NearPixelMinY = 2000;
	dpixel=0;
	dpixelmax=-1;
	pmax=0;
	// Init median calculation Structure. We need the median of values in a 3x3 matrix around current pixel
	Mediator* kmedian = MediatorNew(9);
	
//
//Begin of : Loop all pixels of current frame and search for at least one near pixel
// Remember this is a call back: one pixel per callback is evaluated. A stroke is evaluated in subsequent calls to the callaback
// pixel 0,0		_____________
// pixel 1,0		|0,0|0,1|0,2|
// pixel 2,0		|___|___|___|
// pixel 0,1		|1,0|1,1|1,2|
// pixel 1,1		|___|CUR|___|
// pixel 2,1		|2,0|2,1|2,2|
// pixel 0,2		|___|___|___|
// pixel 1,2
// pixel 2,2

// Initial load of median structure (pixel 1,1 of depth matrix)
// left column
	MediatorInsert(kmedian,t_gamma[depth[0]]);
	MediatorInsert(kmedian,t_gamma[depth[640]]);
	MediatorInsert(kmedian,t_gamma[depth[1280]]);
// center column
	MediatorInsert(kmedian,t_gamma[depth[1]]);
	MediatorInsert(kmedian,t_gamma[depth[641]]);
	MediatorInsert(kmedian,t_gamma[depth[1281]]);
// right column
	MediatorInsert(kmedian,t_gamma[depth[2]]);
	MediatorInsert(kmedian,t_gamma[depth[642]]);
	MediatorInsert(kmedian,t_gamma[depth[1282]]);

	for (i=641; i<FREENECT_FRAME_PIX-640; i++)  // go over inner values od depth matrix (i.e. skip first and last rows and first and last columns)
	{
		if ((i % 640) == 639) i+=2; // if we are at column before last jump to column 1
		//update  median evaluation structure
	// left column
		MediatorInsert(kmedian,t_gamma[depth[i-641]]);
		MediatorInsert(kmedian,t_gamma[depth[i-1]]);
		MediatorInsert(kmedian,t_gamma[depth[i+639]]);
	// center column
		MediatorInsert(kmedian,t_gamma[depth[i-640]]);
		MediatorInsert(kmedian,t_gamma[depth[i]]);
		MediatorInsert(kmedian,t_gamma[depth[i+640]]);
	// right column
		MediatorInsert(kmedian,t_gamma[depth[i-639]]);
		MediatorInsert(kmedian,t_gamma[depth[i+1]]);
		MediatorInsert(kmedian,t_gamma[depth[i+641]]);

	    pval=MediatorMedian(kmedian);
		tx++;
		if(tx >= 640) // end of line reached
		{
			tx = 0;		
			ty++;		// new line
		}
//		Near range pixels : red color

// pval is the median of depth of 9 pixels centered at current pixel
		if (pval < state.config.near_threshold)  // We found a NearPixel
		{
			state.gl_depth_back[3*i+0] = 255; // Red component
			state.gl_depth_back[3*i+1] = 0;   // Green component
			state.gl_depth_back[3*i+2] = 0;   // Blue component
			if (pval > pmax)  // new deeper pixel (i.e. near to screeen)
			{ 
				pmax=pval;
				NearPixelX=tx;
				NearPixelY=ty;
				state.gl_depth_back[3*i+0] = 0;
				state.gl_depth_back[3*i+1] = 255;
			}			
//			if (DistCen[tx][ty]>dpixelmax)  // Current Pixel found is furthest from screen center: update
//			{ 
//				dpixelmax=DistCen[tx][ty];
//				NearPixelX=tx;
//				NearPixelY=ty;
//				state.gl_depth_back[3*i+0] = 0;
//				state.gl_depth_back[3*i+1] = 255;
//			}
			NearPixelCount++;	// we accumulate number of near pixels
			NearPixelXSum+=tx;  // we sum up x coordinates
			NearPixelYSum+=ty;  // we sum up y coordinates
			NearPixelMaxX=MAX(NearPixelMaxX,tx); // we save max x coordinate of nearpixel
			NearPixelMinX=MIN(NearPixelMinX,tx); // we save min x coordinate of nearpixel
			NearPixelMaxY=MAX(NearPixelMaxY,ty); // we save max y coordinate of nearpixel
			NearPixelMinY=MIN(NearPixelMinY,ty); // we save min y coordinate of nearpixel
		}
	
//		Default: mid range pixels : white color
		if (pval >= state.config.near_threshold  && pval < state.config.far_threshold)
		{
			state.gl_depth_back[3*i+0] = 255;
			state.gl_depth_back[3*i+1] = 255;
			state.gl_depth_back[3*i+2] = 255;
		}	
//		Default: far range pixels:  black color
		if (pval >= state.config.far_threshold )
		{
			state.gl_depth_back[3*i+0] = 0;
			state.gl_depth_back[3*i+1] = 0;
			state.gl_depth_back[3*i+2] = 0;
		}
	}

	if(state.config.debug)	printf("Frame Analyzed: NearPixelCount %d\n",NearPixelCount);
//
//End of : Loop all pixels of current frame and search for at least one near pixel
//

//Current Frame evaluated: lets evaluate NearPixels found

// Number of Pixels in Near Blobs is over the threshold NearPixel_TooClose given in input
// This means the subject is too close in current frame
// Reset Pixel Count and restart swipe evaluation
	if(NearPixelCount > state.config.NearPixel_TooClose)   
	{	 
    	if(state.config.debug)	printf("Subject too close\n");
		if(state.config.jsonout && state.config.MMM_Output_status)	printf("{ \"status\" : \"tooclose\"\n}"	);
		state.StrokeEval=1;
	}

// Number of Pixels in Near Blobs is less then threshold NearPixel_TooFarOrNoise given in input but non zero
// This means the subject is within reach but still too far
// Reset Pixel Count and restart swipe evaluation 
	if(NearPixelCount > 0 && NearPixelCount < state.config.NearPixel_TooFarOrNoise)  
	{
    	if(state.config.debug)	printf("Some pixels detected but subject too far\n");
		if(state.config.jsonout && state.config.MMM_Output_status) printf("{ \"status\" : \"somepixels\" }\n" );
		state.StrokeEval=1;
	}

// No near Pixel found
// This means the subject is out of range in current frame
// Reset Pixel Count and restart swipe evaluation 
	if(NearPixelCount == 0)  
	{
    	if(state.config.debug)	printf("Subject too far - Out of reach\n");
		if(state.config.jsonout && state.config.MMM_Output_status)	printf("{ \"status\" : \"toofar\"}\n" );
		state.StrokeEval=1;
	}

// Number of NearPixels in blobs found is neither to small nor too big: subject hand in range
// a swipe is evaluated by evaluating subsequent pixels found between empty frames
// that is : an empty frame (a frame with subject either too far or too close) is considered as a "break" between gestures
// We record x and y coordinates of pixels found in an array
// Begin of section: analyze blob of nearpixels
	if (NearPixelCount !=0 && NearPixelCount < state.config.NearPixel_TooClose && NearPixelCount > state.config.NearPixel_TooFarOrNoise)  
	{		
		state.pointerx = ((NearPixelX-640.0f) / -1); 		// get current x coordinates
		state.pointery = (NearPixelY);					// get current y coordinates
		state.mousex = ((state.pointerx / 630.0f) * state.screenw);	// scale x coordinates to screen size
		state.mousey = ((state.pointery / 470.0f) * state.screenh);	// scale y coordinates to screen size
		mx = state.mousex;			
		my = state.mousey;
    	if(state.config.debug)
    	{ 
    		printf("Subject within range\n");
    		printf("Mouse coordinates  %3d %4d\n",mx,my);
    		printf("Stroke Index  %d \n",state.current_pixel);
    	}
		state.stroke_x[state.current_pixel]=mx; //save current x coord in array of stroke pints
		state.stroke_y[state.current_pixel]=my; //save cuurent y coord in array of stroke pints
		if(state.config.jsonout && state.config.MMM_Output_clicks)	printf("{ \"coord\" : { \"xy\" : \"[ %d , %d]\" }}\n",mx,my );
// This is the first point between invalid or empty frames: reset swipe evaluation support variables
		if(state.current_pixel==0) 
		{  
	    	if(state.config.debug)	printf("First Point in stroke - reset swipe variables\n");
			state.h_stroke_left2right_count=0;
			state.h_stroke_right2left_count=0;
			state.v_stroke_up2down_count=0;
			state.v_stroke_down2up_count=0;
			state.h_stroke_left2right_sum=0;
			state.h_stroke_right2left_sum=0;
			state.v_stroke_up2down_sum=0;
			state.v_stroke_down2up_sum=0;
			state.h_sum=0;
			state.v_sum=0;
			state.h_mean=mx;
			state.v_mean=my;
			state.h_variance=0;
			state.v_variance=0;
			state.left2rightmul=0;
			state.right2leftmul=0;
			state.up2downmul=0;
			state.down2upmul=0;
			state.StrokeEval=0;
		} 
		else // this is an additional point in a stroke. Update swipe statistics variables
		{  
    		if(state.config.debug)	printf("Update  swipe variables\n");
    		if (state.stroke_x[state.current_pixel]>state.stroke_x[state.current_pixel-1]) 
    		{ 
    			state.h_stroke_left2right_count++;
    			state.h_stroke_left2right_sum += abs(state.stroke_x[state.current_pixel]-state.stroke_x[state.current_pixel-1]);
    		}  
    		else 
    		{
    			state.h_stroke_right2left_count++;
    			state.h_stroke_right2left_sum += abs(state.stroke_x[state.current_pixel-1]-state.stroke_x[state.current_pixel]);
    		}
    		if(state.stroke_y[state.current_pixel]>state.stroke_y[state.current_pixel-1]) 
    		{
    			state.v_stroke_up2down_count++;
    			state.v_stroke_up2down_sum+=abs(state.stroke_y[state.current_pixel-1]-state.stroke_y[state.current_pixel]);
    		}  
    		else 
    		{	
    			state.v_stroke_down2up_count++;
    			state.v_stroke_down2up_sum+=abs(state.stroke_y[state.current_pixel]-state.stroke_y[state.current_pixel-1]);
    		}
		}
    	state.h_sum+=state.stroke_x[state.current_pixel];
   		state.v_sum+=state.stroke_y[state.current_pixel];
   		if(state.config.debug)
   		{ 
   			printf("Deltas: \tH %4d \tV %4d\n",abs(state.stroke_x[state.current_pixel-1]-state.stroke_x[state.current_pixel]),abs(state.stroke_y[state.current_pixel-1]-state.stroke_y[state.current_pixel]));
   			printf("Counts: \tL2R %5d\tR2L %5d\tU2D %5d\tD2U %5d\n",state.h_stroke_left2right_count,state.h_stroke_right2left_count,state.v_stroke_up2down_count,state.v_stroke_down2up_count);
   			printf("StrokeSums: \tL2R %5d\tR2L %5d\tU2D %5d\tD2U %5d\n",state.h_stroke_left2right_sum,state.h_stroke_right2left_sum,state.v_stroke_up2down_sum,state.v_stroke_down2up_sum);
   			printf("HSUm VSUM: \tHsum %10d\tVSum %10d\n",state.h_sum,state.v_sum);
   		}

// If current evaluated pixel coordinates are within square area defined by input parameter gesture_click_area
// Increment the hovering counter and reset Stroke index  		StrokeSums
		if ((state.PointerX <= (mx + state.config.gesture_click_area))  && (state.PointerX >= (mx -state.config.gesture_click_area)) && (state.PointerY <= (my + state.config.gesture_click_area))  && (state.PointerY >= (my - state.config.gesture_click_area))) 
		{
			state.current_hovering_cycles++;
	    	if(state.config.debug)	printf("Mouse Hovering : Current hovering cycle %3d\n",state.current_hovering_cycles);
		} 
		else  
// Current evaluated pixel coordinates are not within square area defined by input parameter gesture_click_area
// Reset the hovering counter and increment stroke pixel index   			
		{
			state.PointerX = mx; //New initial position X
			state.PointerY = my; //New initial position Y
			state.current_hovering_cycles = 0; // Restart counting current_hovering_cycles
		}		
		state.current_pixel++;
// Check if mouse was hovering for more then hovering_threshold subsequent frames over the click area
// Simulate click at the point
// Set debounce count to avoid double clicks    			
		if(state.current_hovering_cycles > state.config.hovering_threshold) 
		{
			state.current_hovering_cycles = -state.config.hovering_threshold*2;  		// set debounce count
			XTestFakeButtonEvent(state.display, 1, TRUE, CurrentTime);  	// send mouse lmb down 
			XTestFakeButtonEvent(state.display, 1, FALSE, CurrentTime);	// send mouse lmb up
			state.current_pixel=0;  //Reset Stroke Pixel Index 
			state.StrokeEval=0;
	    	if(state.config.debug)	printf("Click at \tX %d\tY %dn",mx,my);
			if(state.config.jsonout && state.config.MMM_Output_clicks)	printf("{ \"click\" : { \"xy\" : \"[ %d , %d]\" }}\n",mx,my );
		}
		XTestFakeMotionEvent(state.display, -1, mx, my, CurrentTime);		// send mouse movement
		if(state.config.debug)	printf("Coordinates \tX %3d\tY %3d\n",mx,my);
//		if(state.jsonout && state.MMM_Output_coords)	printf("{ \"coords\" : { \"xy\" : \"[ %d , %d]\" }}\n",mx,my );
		XSync(state.display, 0);
	}
	// End of section: analyze blob of nearpixels

// begin of section: Evaluate Swipe if frame empty or not in threshold 
	if(state.StrokeEval)
	{
// Begin of section: We have a sequence of points that are between the allowed paramenters set
		if (state.current_pixel>state.config.minimum_stroke_points && state.current_pixel<state.config.maximum_stroke_points) 
		{
			state.h_mean=state.h_sum/state.current_pixel;  //calculate statistics for euristic
			state.v_mean=state.v_sum/state.current_pixel;
			for (j=0;j<state.current_pixel;j++)	
			{
				state.h_variance+=pow(state.stroke_x[j]-state.h_mean,2);
				state.v_variance+=pow(state.stroke_y[j]-state.v_mean,2);
				if(state.config.debug) printf("N X Y \t%3d\t%3d\t%4d\n",j, state.stroke_x[j],state.stroke_y[j]);
			}
			state.h_variance=sqrt(state.h_variance/state.current_pixel);
			state.v_variance=sqrt(state.v_variance/state.current_pixel);
			state.left2rightmul = state.h_stroke_left2right_count*state.h_stroke_left2right_sum;
			state.right2leftmul = state.h_stroke_right2left_count*state.h_stroke_right2left_sum;
			state.up2downmul = state.v_stroke_up2down_count*state.v_stroke_up2down_sum;
			state.down2upmul = state.v_stroke_down2up_count*state.v_stroke_down2up_sum;
			if(state.config.debug)
			{ 
				printf("____________________________________________\n");
				printf("state.current_pixel, %i\n",state.current_pixel);
				printf("left2rightcnt\tright2leftcnt\t%8i\t%8i\n",state.h_stroke_left2right_count, state.h_stroke_right2left_count);
				printf("up2downcnt\tdown2upcnt\t%8i\t%8i\n\n",state.v_stroke_up2down_count,state.v_stroke_down2up_count);
				printf("left2rightsum\tright2leftsum\t%8i\t%8i\n",state.h_stroke_left2right_sum, state.h_stroke_right2left_sum);
				printf("up2downsum\tdown2up_sum\t%8i\t%8i\n\n",state.v_stroke_up2down_sum,state.v_stroke_down2up_sum);
				printf("left2rightmul\tright2leftmul\t%20d\t%20d\n",state.left2rightmul,state.right2leftmul);
				printf("up2downmul\tdown2upmul\t%20d\t%20d\n\n",state.up2downmul,state.down2upmul);
				printf("h_mean\tv_mean\t%5.5f\t%5.5f\n",state.h_mean,state.v_mean);
				printf("h_variance\tv_variance\t%5.5f\t%5.5f\n",state.h_variance,state.v_variance);
				printf("____________________________________________\n");
			}
			if (state.h_variance<state.config.h_varmax || state.v_variance<state.config.v_varmax)  // Either Variance is belo variance threshold parameter
			{
				if (state.h_variance<state.config.h_varmax && state.v_variance<state.config.v_varmax)
				{ 
					if(state.config.debug) printf("Garbage ");
				}
				else
				{ 
					if(state.h_variance<state.v_variance)
				   	{ 
						if(state.config.debug) if(state.up2downmul>state.down2upmul) printf("Down \n"); else printf("Up \n");				
						if(state.config.jsonout && state.config.MMM_Output_swipes) if(state.up2downmul>state.down2upmul) printf("{ \"swipe\" : \"down\" }\n" ); else printf("{ \"swipe\" : \"up\" }\n" );				
				   	}
				   	else
				   	{ 
				   		if(state.config.debug) if(state.left2rightmul>state.right2leftmul) printf("right \n"); else printf("left \n");				
						if(state.config.jsonout && state.config.MMM_Output_swipes) if(state.left2rightmul>state.right2leftmul)	printf("{ \"swipe\" : \"right\" }\n" ); else printf("{ \"swipe\" : \"left\" }\n" );				
				   	}
				}
			}	
			else
			{
					if(state.config.debug) printf("No Swipe\n");
			}	
		}  // End of section: We have a sequence of points that are between the allowed paramenters 
		if(state.config.debugstop) getchar();
		state.StrokeEval=0;
		state.current_pixel=0;
	}
	// end of section: Evaluate Swipe if frame empty or not in threshold 
		
	state.got_frames++;
	pthread_cond_signal(&state.gl_frame_cond);
	pthread_mutex_unlock(&state.gl_backbuf_mutex);
	if(state.config.debug) printf("___________________________ENDOFRAME_________________________\n\n");
}

void rgb_cb(freenect_device *dev, void *rgb, uint32_t timestamp)
{
	pthread_mutex_lock(&state.gl_backbuf_mutex);
	state.got_frames++;
	memcpy(state.gl_rgb_back, rgb, FREENECT_VIDEO_RGB_SIZE);
	pthread_cond_signal(&state.gl_frame_cond);
	pthread_mutex_unlock(&state.gl_backbuf_mutex);
}

void *freenect_threadfunc(void *arg)
{
	freenect_set_tilt_degs(state.f_dev, state.config.freenect_angle);
	freenect_set_led(state.f_dev, state.config.freenect_led);
	freenect_set_depth_callback(state.f_dev, depth_cb);
	freenect_set_video_callback(state.f_dev, rgb_cb);
	freenect_set_video_buffer(state.f_dev, FREENECT_VIDEO_RGB);  
	freenect_set_depth_buffer(state.f_dev, FREENECT_DEPTH_11BIT);

	freenect_start_depth(state.f_dev);
	freenect_start_video(state.f_dev);

	//printf("'W'-Tilt Up, 'S'-Level, 'X'-Tilt Down, '0'-'6'-LED Mode\n");

	while(!state.die && freenect_process_events(state.f_ctx) >= 0 )
	{
		freenect_raw_tilt_state* tilt_state;
		freenect_update_tilt_state(state.f_dev);
		tilt_state = freenect_get_tilt_state(state.f_dev);;
		double dx,dy,dz;
		freenect_get_mks_accel(tilt_state, &dx, &dy, &dz);
		//printf("\r raw acceleration: %4d %4d %4d  mks acceleration: %4f %4f %4f\r", ax, ay, az, dx, dy, dz);
		fflush(stdout);
	}

	
	if(state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"Start Shutting Down Streams\"}\n");
	if(state.config.debug) printf("Start Shutting Down Streams");

	freenect_stop_depth(state.f_dev);
	freenect_stop_video(state.f_dev);

	freenect_close_device(state.f_dev);
	freenect_shutdown(state.f_ctx);
	if(state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"Done Shutting Down Streams\"}\n");
	if(state.config.debug) printf("Done Shutting Down Streams");
	return NULL;
}


int main(int argc, char **argv)
{
	int res;

    if ((argc != 25) || ((argc == 1) && strcmp (argv[1],"--help")))
	{
		if (argc!=25)
			if (state.config.jsonout && state.config.MMM_Output_log) 	printf("{ \"log\" : \"Wrong Number of Parameters: %2d \"}\n",argc);
		printf("Number of Parameters %2d \n",argc);
		printf("- NearPixel_TooClose: Number of maximum near pixel to accept before sending a too close message\n");
		printf("- NearPixel_TooFarOrNoise: Number of minimum near pixel to accept as pointer (i.e. not noise)\n");
		printf("- KinectLogLevel: Log level to output (0-7) 0 = nothing / 0 Flood\n");
		printf("- KinectSwhowScreen: 0: no output 1: Show camera and depth camera\n");
		printf("- KinectAngle: Start angle for kinect -30 : 30\n");
		printf("- KinectLed: Led color to light on 0-6\n");
		printf("- ClickSize: Size of the area to be considered as mouse steady for click\n");
		printf("- ClickPauseCount: number of subsequent frames with steady mouse to trigger a click \n");
		printf("- minimum_stroke_points: minimum number of points to evaluate as a stroke\n");
		printf("- maximum_stroke_points: maximum number of points to evaluate as a stroke\n");
		printf("- Horizontal Variance threshold: maximum variance of coords in horizontal direction for a set of coords to be considered a vertical swipe\n");
		printf("- Vertical Variance threshold: maximum variance of coords in vertical direction for a set of coords to be considered an horizontal swipe\n");
		printf("- near_threshold: depth for near points\n");
		printf("- far_threshold: depth for far points\n");
		printf("- Elbow X\n");
		printf("- ElbowY\n");
		printf("- JSon Output X\n");
		printf("- Output program log\n");
		printf("- Output sensor status\n");
		printf("- Output Click events\n");
		printf("- Output Coordinates events\n");
		printf("- Output Swipe events\n");
		printf("- Verbose debug to stdout\n");
		printf("- Pause Verbose debug output at swipe evaluation\n");
		return 1;
	}
	else	{
//			printf("{ \"info\" : \"NearPixel_TooClose: Number of maximum near pixel to accept before sending a too close message %i\n",NearPixel_TooClose\"}\n");
//			printf("{ \"info\" : \"NearPixel_TooFarOrNoise: Number of minimum near pixel to accept as pointer (i.e. not noise) %i\n",NearPixel_TooFarOrNoise\"}\n");
//			printf("{ \"info\" : \"KinectLogLevel: Log level to output (0-7) 0 = nothing / 0 Flood %i\n",freenect_log_level\"}\n");
//			printf("{ \"info\" : \"KinectSwhowScreen: 0: no output 1: Show camera and depth camera %i\n",debug\"}\n");
//			printf("{ \"info\" : \"KinectAngle: Start angle for kinect -30 : 30  %i\n",freenect_angle\"}\n");
//			printf("{ \"info\" : \"KinectLed: Led color to light on 0-6  %i\n",freenect_led\"}\n");
//			printf("{ \"info\" : \"ClickSize: Size of the area to be considered as mouse steady for click %i\n",gesture_click_area\"}\n");
//			printf("{ \"info\" : \"ClickPauseCount: number of subsequent frames with steady mouse to trigger a click %i\n",hovering_threshold\"}\n");
//			printf("{ \"info\" : \"minimum_stroke_points: minimum number of points to evaluate as a stroke %i\n",minimum_stroke_points\"}\n");
//			printf("{ \"info\" : \"maximum_stroke_points: maximum number of points to evaluate as a stroke %i\n",maximum_stroke_points\"}\n");
//			printf("{ \"info\" : \"near_threshold: depth for near points %i\n",near_threshold\"}\n");
//			printf("{ \"info\" : \"far_threshold: depth for far points %i\n",far_threshold\"}\n");
//			printf("{ \"info\" : \"ystretch %f\n\n",ystretch\"}\n");
//			printf("{ \"info\" : \"	RorriM
//			printf("{ \"info\" : \"ElbowY %i\n", ScreenCenterX\"}\n");
//			printf("{ \"info\" : \"ElbowY %i\n", ScreenCenterY\"}\n");
//			printf("{ \"info\" : \"JSon Output %d\n", jsonout\"}\n");
//			printf("{ \"info\" : \"Send Click events %d\n", MMM_Output_clicks\"}\n");
//			printf("{ \"info\" : \"Send Coordinates events %d\n", MMM_Output_coords\"}\n");
//			printf("{ \"info\" : \"Send Swipe events %d\n", MMM_Output_swipes\"}\n");
//			printf("{ \"info\" : \"Verbose debug %d\n", debug\"}\n");

		state.config.NearPixel_TooClose = atoi(argv[1]);
		state.config.NearPixel_TooFarOrNoise = atoi(argv[2]);
		state.config.freenect_log_level = atoi(argv[3]); 
		state.config.ShowScreen = atoi(argv[4]);
		state.config.freenect_angle = atoi(argv[5]);
		state.config.freenect_led = atoi(argv[6]);
		state.config.gesture_click_area =atoi(argv[7]);
		state.config.hovering_threshold = atoi(argv[8]);
		state.config.minimum_stroke_points = atoi(argv[9]);
		state.config.maximum_stroke_points = atoi(argv[10]);
		state.config.h_varmax = atoi(argv[11]);
		state.config.v_varmax = atoi(argv[12]);
		state.config.near_threshold = atoi(argv[13]);
		state.config.far_threshold = atoi(argv[14]);
		state.config.ScreenCenterX = atoi(argv[15]);
		state.config.ScreenCenterY = atoi(argv[16]); 
		state.config.jsonout = atoi(argv[17]);
		state.config.MMM_Output_status =  atoi(argv[18]);
		state.config.MMM_Output_log = atoi(argv[19]);
		state.config.MMM_Output_clicks =  atoi(argv[20]);
		state.config.MMM_Output_coords = atoi(argv[21]);
		state.config.MMM_Output_swipes = atoi(argv[22]); 
		state.config.debug = atoi(argv[23]); 
		state.config.debugstop = atoi(argv[24]); 

		if((state.config.jsonout && state.config.MMM_Output_log) || state.config.debug)	{
			printf("{ \"log\" : \"Kinect Mouse and Swipe starting\"}\n");
			printf("{ \"log\" : \"Parsing Args\"}\n");
			printf("{ \"log\" : \"NearPixel_TooClose %i \"}\n",state.config.NearPixel_TooClose);
			printf("{ \"log\" : \"NearPixel_TooFarOrNoise %i \"}\n",state.config.NearPixel_TooFarOrNoise);
			printf("{ \"log\" : \"KinectLogLevel %i \"}\n",state.config.freenect_log_level);
			printf("{ \"log\" : \"KinectSwhowScreen %i \"}\n",state.config.ShowScreen);
			printf("{ \"log\" : \"KinectAngle %i \" }\n",state.config.freenect_angle);
			printf("{ \"log\" : \"KinectLed %i \" }\n",state.config.freenect_led);
			printf("{ \"log\" : \"ClickSize %i \" }\n",state.config.gesture_click_area);
			printf("{ \"log\" : \"ClickPauseCount %i \"}\n",state.config.hovering_threshold);
			printf("{ \"log\" : \"minimum_stroke_points %i \"}\n",state.config.minimum_stroke_points);
			printf("{ \"log\" : \"maximum_stroke_points %i \"}\n",state.config.maximum_stroke_points);
			printf("{ \"log\" : \"horizontal variance threshold for swipe %i \"}\n",state.config.minimum_stroke_points);
			printf("{ \"log\" : \"vertical variance threshold for swipe%i \"}\n",state.config.maximum_stroke_points);
			printf("{ \"log\" : \"near_threshold: %i \"}\n" ,state.config.near_threshold);
			printf("{ \"log\" : \"far_threshold: depth for far points %i \"}\n",state.config.far_threshold);
			printf("{ \"log\" : \"ElbowX %i \"}\n",state.config.ScreenCenterX);
			printf("{ \"log\" : \"ElbowY %i \"}\n", state.config.ScreenCenterY);
			printf("{ \"log\" : \"JSon Output %d \"}\n", state.config.jsonout);
			printf("{ \"log\" : \"Output status to stdout %d \"}\n", state.config.MMM_Output_status);
			printf("{ \"log\" : \"Output log to stdout %d \"}\n", state.config.MMM_Output_log);
			printf("{ \"log\" : \"Output Click events to stdout %d \"}\n", state.config.MMM_Output_clicks);
			printf("{ \"log\" : \"Output Coordinates events to stdout %d \"}\n", state.config.MMM_Output_coords);
			printf("{ \"log\" : \"Output Swipe events to stdout %d \"}\n", state.config.MMM_Output_swipes);
			printf("{ \"log\" : \"Verbose debug %d \"}\n", state.config.debug);
			printf("{ \"log\" : \"Verbose debug stop at swipe eval %d \"}\n", state.config.debugstop);
		}
	}
	
	
	//mousemask(ALL_MOUSE_EVENTS, NULL);
	if(state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"Opening Display\"}\n");
	if(state.config.debug) printf("Opening display \n");

	state.display = XOpenDisplay(0);

	// Careful: this will dump if not run from terminal directly on system running XServer 
	state.root_window = (state.display);

	state.screenw = XDisplayWidth(state.display, SCREEN);
	state.screenh = XDisplayHeight(state.display, SCREEN);

	if(state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"Default Display Found\"}\n{ \"log\" : \"Display Size %d %d }\n", state.screenw, state.screenh);
	if(state.config.debug) printf("Default Display Found\nDisplay Size %d %d \n", state.screenw, state.screenh);

//	state.screenw += 200;
//	state.screenh += 200;

	int i,j;
	for (i=0; i<2048; i++) {
		float v = i/2048.0;
		v = powf(v, 3)* 6;
		t_gamma[i] = v*6*256;
	}

	state.g_argc = argc;
	state.g_argv = argv;

	/* Precalculate all distances from Screen Center */
	for (i=0; i<640;i++)
			for (j=0; j<480;j++)
				state.DistCen[i][j] = sqrt(pow(i-state.config.ScreenCenterX,2)+pow(j-state.config.ScreenCenterY,2));

	
	if (freenect_init(&state.f_ctx, NULL) < 0) {
		if (state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"freenect_init() failed\" }\n");
		if (state.config.debug) printf("Error freenect_init() failed\n");
		return 1;
	}

	freenect_set_log_level(state.f_ctx, state.config.freenect_log_level);

	int nr_devices = freenect_num_devices (state.f_ctx);
		if (state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"%d devices Found\" }\n", nr_devices);
		if (state.config.debug) printf("%d devices Found\n", nr_devices);
	

	int user_device_number = 0;

	if (freenect_open_device(state.f_ctx, &state.f_dev, user_device_number) < 0) {
		if (state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"No Kinect found\" \n}");
		if (state.config.debug) printf("Error : No Kinect found\n");
		return 1;
	}

	res = pthread_create(&state.freenect_thread, NULL, freenect_threadfunc, NULL);
	if (res) {
		if (state.config.jsonout && state.config.MMM_Output_log) printf("{ \"log\" : \"Could not create thread\" }\n");
		if (state.config.debug) printf("Error could not create thread.\n");
		return 1;
	}

	gl_threadfunc(NULL);

	return 0;
}
