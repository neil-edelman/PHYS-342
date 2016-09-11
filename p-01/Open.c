#include <stdlib.h> /* malloc free */
#include <stdio.h>  /* fprintf */
#ifdef GL
#include <GL/gl.h>
#include <GL/glut.h>
#else
#include <OpenGL/gl.h> /* OpenGL ** may be GL/gl.h */
#include <GLUT/glut.h> /* GLUT */
#endif
#include "Simulation.h"

struct Open *Open(const int width, const int height, const char *title);
void Open_(void);
void update(int);
void display(void);
void resize(int width, int height);
void keyUp(unsigned char k, int x, int y);
void keyDn(unsigned char k, int x, int y);
void keySpecial(int key, int x, int y);
void cube(const float x, const float y, const float z, const float a);

struct Open {
	struct Simulation *s;
	float             rot;
	int               frame;
};

struct Open *open = 0;

const static int granularity   = 32;
const static float speed       = .7;
static float black_of_space[4] = { 0, 0, 0, 0 }/*{ 1, 1, 1, 0 }*/;
static float current[4]        = { 0, 0, 1, 1 }/*{ .4, .4, .7, 1 }*/;
static float magnetic[4]       = { .2, 1, .9, .3 }/*{ .2, .2, .3, .3 }*/;

int main(int argc, char **argv) {
	/* negotiate with library */
	glutInit(&argc, argv);

	if(!Open(360, 240, "Simuation")) return EXIT_FAILURE;
	/* atexit because the loop never returns */
	if(atexit(&Open_)) perror("~Open");
	glutMainLoop();

	return EXIT_SUCCESS;
}

struct Open *Open(const int width, const int height, const char *title) {
	GLfloat lightPos[4] = { 1.0, 10.0, 10.0, 0.0 }, lightAmb[4] = { 1, .5, .2, 1 };

	if(open || width <= 0 || height <= 0 || !title) {
		fprintf(stderr, "Open: error initialising.\n");
		return 0;
	}
	if(!(open = malloc(sizeof(struct Open)))) {
		perror("Open constructor");
		Open_();
		return 0;
	}
	open->s     = 0;
	open->rot   = 0;
	open->frame = 0;	
	if(!(open->s = Simulation(granularity, &glVertex3f))) { Open_(); return 0; }
	fprintf(stderr, "Open: new, #%p.\n", (void *)open);
	/* initial conditions */
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH); /* RGB[A] is implied */
	glutInitWindowSize(width, height); /* just a suggestion */
	/* create */
	glutCreateWindow(title);
	/* initialise */
	glShadeModel(GL_SMOOTH);
	glClearColor(black_of_space[0], black_of_space[1], black_of_space[2], black_of_space[3]);
	glClearDepth(1.0);
	/*glEnable*/glDisable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	/*glEnable(GL_LIGHTING);*/
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (GLfloat *)&lightAmb);
	glShadeModel(GL_FLAT);
	glutReshapeWindow(width, height);
	/* set callbacks */
	glutDisplayFunc(&display);
	glutReshapeFunc(&resize);
	glutKeyboardFunc(&keyDn);
	/* glutIdleFunc(0); disable */
	glutTimerFunc(25, update, 0);

	return open;
}

void Open_(void) {
	if(!open) return;
	fprintf(stderr, "~Open: erase, #%p.\n", (void *)open);
	if(open->s) Simulation_(&open->s);
	free(open);
	open = 0;
}

/* private */

void update(int value) {
	int (*e)(struct Simulation *, int);

	open->rot += speed;
	glutPostRedisplay();
	glutTimerFunc(25, update, 0);
	SimulationUpdate(open->s);
	if((e = SimulationGetExplode(open->s))) {
		if((e(open->s, open->frame))) {
			open->frame++;
		} else {
			SimulationClearExplode(open->s);
			open->frame = 0;
		}
	}
}

void display(void) {
	int size;
	GLfloat x, y, z, a, offset;

	/* clear screen and depthbuf, make sure it's modelview, and reset the matrix */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	size   = SimulationGetSize(open->s);
	offset = -(float)size / 2 + .5;
	glTranslatef(0, 0, 3 * offset);
	glRotatef(open->rot, 0, 1, 0);
	glTranslatef(offset, offset + .2, offset);
	glBegin(GL_LINES);
	glColor4f(current[0], current[1], current[2], current[3]);
	SimulationCurrent(open->s);
	glColor4f(magnetic[0], magnetic[1], magnetic[2], magnetic[3]);
	SimulationMagnetic(open->s);
	glEnd();
	for(x = 0; x < size; x += 1) {
		for(y = 0; y < size; y += 1) {
			for(z = 0; z < size; z += 1) {
				a = 1 - SimulationGetMu(open->s, x, y, z);
				glPointSize(a * 64 + 1);
				/*glColor4f(.5*(1 - a), .2, .5*(1 - a), .3);*/
				glColor4f(1. - a, 5., 1. - a, .3);
				glBegin(GL_POINTS);
				glVertex3f(x + .5, y + .5, z + .5);
				glEnd();
			}
		}
	}
	glutSwapBuffers();
}

void resize(int width, int height) {
	if(width <= 0 || height <= 0) return;
	glViewport(0, 0, width, height);
	/* calculate the projection */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f,(float)width / height, 0.1f, 100.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, -500/*-23.0 / 2.5*/);
}

void keyDn(unsigned char k, int x, int y) {
	SimulationAnimation(open->s);
}
