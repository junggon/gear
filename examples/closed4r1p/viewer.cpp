/* modified simpleviewer.cpp (http://www.gvu.gatech.edu/~jarek/courses/openGL/viewer.htm) */
/* Simple geometry viewer:  Left mouse: rotate;  Middle mouse:  zoom;  Right mouse:   menu;  ESC to quit */
#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <GL/freeglut.h>
#else
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif
#endif
#include "myworld.h"

#define KEY_SPACE 32
#define KEY_ESC 27
#define KEY_ENTER 13


// mysystem
MyWorld *_pmyworld = NULL;
float _mscale = 1; // model scale 

// simulation
double _step_size=0.001;
double _current_time=0.0; // this is also updated during replay
double _max_simul_time=5;
int _simul_save_freq = 120; // data saving frequency
int _max_counter_save = int(1./(double(_simul_save_freq)*_step_size));

// replay simulated motion
bool _breplay=false; 
int _current_frame_idx=0; // current replay frame index
double _speed=1.0; // replay _speed
int _replay_render_freq=33; // rendering frequency when replay
int _num_skip_frames = int(double(_simul_save_freq)/double(_replay_render_freq)*_speed-1); // number of frames to be skipped when replay

// system color
GLfloat _textcolor[] = {0.0f,0.0f,0.0f,1.0f};
GLfloat _backgroundcolor[] = {1.0f,1.0f,1.0f,1.0f};

// viewer state
Vec3 _spoint_prev;
SE3 _mvmatrix, _mvmatrix_prev; // view change
double &_xcam = _mvmatrix[12], &ycam = _mvmatrix[13]; // view translate
double &_sdepth = _mvmatrix[14]; // zoom in/out
float _fovy=64.0f, _zNear=0.001f, _zFar=10000.0f;
float _aspect = 5.0f/4.0f;
long _winsizeX, _winsizeY;
int _downX, _downY;
bool _bRotateView = false, _bTranslateView = false, _bZoomInOut = false, _bSliderControl = false, _bShowUsage = false;
GLfloat _light0Position[] = {0,1,0,1.0};

using namespace std;

void TestFunc(void)
{
	_pmyworld->testfunc();
}

//ofstream fout("energy.m");
void RunSimulation(void)
{
	int cnt = 0;
	MyWorld::SimulData sd;
	_current_time = 0;
	_pmyworld->simuldata.clear();

	// save initial state
	sd.t = _current_time;
	sd.q = _pmyworld->system.get_q();
	sd.dq = _pmyworld->system.get_dq();
	_pmyworld->simuldata.push_back(sd);
	cnt = 0;

	tic();
	while ( _current_time <= _max_simul_time ) {
		// init body force and joint  torque
		_pmyworld->system.initBodyForcesAndJointTorques();

		// run forward dynamics simulation
		if ( !_pmyworld->system.stepSimulation(_step_size) ) {
			cout << "error:: failed in stepping simulation forward at t = " << _current_time << "sec" << endl;
			break;
		}

		// increase time
		_current_time += _step_size;

		// save simulation data
		if ( cnt++ >= _max_counter_save ) {
			sd.t = _current_time;
			sd.q = _pmyworld->system.get_q();
			sd.dq = _pmyworld->system.get_dq();
			_pmyworld->simuldata.push_back(sd);
			cnt = 0;
			//double e1 = _pmyworld->system.getGravitationalPotentialEnergy();
			//double e2 = _pmyworld->system.getKineticEnergy();
			//fout << _current_time << " " << e1+e2 << " " << e1 << " " << e2 << endl;
		}
	}
	double t_elapsed = toc();

	cout << "elapsed time for " << _max_simul_time << " sec simulation = " << t_elapsed << " sec" << endl;
}

void ReadSimulData(void)
{
	if ( _pmyworld->simuldata.size() == 0 ) return;

	if ( _current_frame_idx >= _pmyworld->simuldata.size() ) {
		_current_frame_idx = _pmyworld->simuldata.size()-1;
	}
	if ( _current_frame_idx < 0 ) {
		_current_frame_idx = 0;
	}

	_current_time = _pmyworld->simuldata[_current_frame_idx].t;
	_pmyworld->system.set_q(_pmyworld->simuldata[_current_frame_idx].q);
	_pmyworld->system.set_dq(_pmyworld->simuldata[_current_frame_idx].dq);
	_pmyworld->system.updateKinematics();
}

void PlotString(void *font, int x, int y, const char *string)
{        
#ifndef __APPLE__
	// switch to projection mode
	glMatrixMode(GL_PROJECTION);	
	glPushMatrix();	
	// reset matrix	
	glLoadIdentity();	
	//Set the viewport	
	glViewport(0, 0, _winsizeX, _winsizeY);
	// set a 2D orthographic projection	
	glOrtho(0, _winsizeX, _winsizeY, 0, -100, 100);	
	glMatrixMode(GL_MODELVIEW);	
	glPushMatrix();	
	glLoadIdentity();
	glTranslatef(0.0f,0.0f,99.9f);
	// draw
	glRasterPos2i(x,y);
	glutBitmapString( font, (const unsigned char *)string );
	// restore
	glMatrixMode(GL_PROJECTION);	
	glPopMatrix ();	
	glMatrixMode(GL_MODELVIEW);	
	glPopMatrix ();
#endif
}

void DrawTextInfo(void)
{
	glDisable(GL_LIGHTING);
	char *str = new char[150]; 
	sprintf(str, "t = %6.3f", _current_time);
	char str2[] = "ESC : exit, ENTER : start simulation, SPACE : start/pause replay, t: run test function";
	glColor4fv(_textcolor);
	PlotString(GLUT_BITMAP_HELVETICA_18, 20, 20, str);
	PlotString(GLUT_BITMAP_HELVETICA_12, 20, _winsizeY-20, str2);
	glEnable(GL_LIGHTING);
}

void StopReplay(void)
{
	_breplay = false;
}

void StartReplay(void)
{
	_breplay = true;
}

void DrawModel(void) 
{
	// draw coordinate frame
	GLfloat lw;
	glGetFloatv(GL_LINE_WIDTH, &lw);
	glLineWidth (1.0);
	glBegin (GL_LINES);
	glColor3f (1,0,0); // X axis is red.
	glVertex3f(0,0,0);
	glVertex3f(1,0,0);
	glColor3f (0,1,0); // Y axis is green.
	glVertex3f(0,0,0);
	glVertex3f(0,1,0);
	glColor3f (0,0,1); // z axis is blue.
	glVertex3f(0,0,0);
	glVertex3f(0,0,1);
	glEnd();
	glLineWidth(lw);

	// draw mysystem
	_pmyworld->render();
}

void Timer(int extra) 
{
	if ( _breplay ) {
		_current_frame_idx += 1 + _num_skip_frames;
		if ( _current_frame_idx >= (int)_pmyworld->simuldata.size() ) {
			_breplay = false;
			_current_frame_idx = (int)_pmyworld->simuldata.size()-1;
		}
		ReadSimulData();
		glutPostRedisplay();
	}
	glutTimerFunc(_replay_render_freq, Timer, 0);
}

void ResetTime(void) 
{
	_current_frame_idx = 0;
	ReadSimulData();
	glutPostRedisplay();
}

void ReshapeCallback(int width, int height) 
{
	_winsizeX = width;
	_winsizeY = height;
	_aspect = (float)_winsizeX/(float)_winsizeY;
	glViewport(0, 0, _winsizeX, _winsizeY);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glutPostRedisplay();
}
 
void DisplayCallback(void) 
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(_fovy, _aspect, _zNear, _zFar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(_mvmatrix.GetArray());
	DrawModel();
	DrawTextInfo();
	glutSwapBuffers();
	glutPostRedisplay(); 
	glClearColor(_backgroundcolor[0],_backgroundcolor[1],_backgroundcolor[2],_backgroundcolor[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
}

void KeyboardCallback(unsigned char ch, int x, int y) 
{
	switch (ch) {
		case KEY_ESC:// exit program
			std::cout << "exit program" << std::endl;
			exit(0);
			break;
		case KEY_ENTER:  // start/pause simulation
			RunSimulation();
			_current_frame_idx = _pmyworld->simuldata.size()-1;
			ReadSimulData();
			break;
		case KEY_SPACE: // play/pause replay
			if ( _breplay ) {
				StopReplay();
			} else {
				if ( _current_frame_idx == _pmyworld->simuldata.size()-1 ) {
					ResetTime();
				}
				StartReplay();
			}
			break;
		case 'r': // reset time
			ResetTime();
			break;
		case 't': // test function
			TestFunc();
			break;
	}
	glutPostRedisplay(); 
}

void SpecialKeyboardCallback(int key, int x, int y) 
{
	switch (key) {
		case GLUT_KEY_RIGHT:
			break;
		case GLUT_KEY_LEFT:
			break;
		case GLUT_KEY_DOWN:
			break;
		case GLUT_KEY_UP:
			break;
	} 
	glutPostRedisplay(); 
}
 
Vec3 get_pos_sp(int x, int y)
{
	Vec3 p;
	p[0] = 2.0 * double(x) / double(_winsizeX) - 1.0;
	p[1] = 1.0 - 2.0 * double(y) / double(_winsizeY);
	p[2] = p[0] * p[0] + p[1] * p[1];
	if ( p[2] < 1.0 ) p[2] = sqrt(1.0 - p[2]);
	else { p[2] = sqrt(p[2]); p[0] /= p[2]; p[1] /= p[2]; p[2] = 0.0; }
	return p;
}

void MouseCallback(int button, int state, int x, int y) 
{
	_downX = x; _downY = y;
	_bRotateView = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN) && y < _winsizeY-50 );
	_bTranslateView = ((button == GLUT_MIDDLE_BUTTON) &&  (state == GLUT_DOWN));
	_bZoomInOut = ((button == GLUT_RIGHT_BUTTON) &&  (state == GLUT_DOWN));

	if ( _bRotateView ) {
		_mvmatrix_prev = _mvmatrix;
		_spoint_prev = get_pos_sp(x, y);
	}

	// wheel (zoom in/out)
	if ( (button != GLUT_LEFT_BUTTON) && (button != GLUT_MIDDLE_BUTTON) && (button != GLUT_RIGHT_BUTTON) ) {
		if ( button==3 ) { // wheel up --> zoom out
			_sdepth -= 0.1f*_mscale;
		} 
		if ( button==4 ) { // wheel down --> zoom in
			_sdepth += 0.1f*_mscale;
		}
	}

	glutPostRedisplay();
}
 
void MotionCallback(int x, int y) 
{
	bool bctrlpressed = false;
	if (glutGetModifiers() == GLUT_ACTIVE_CTRL) { 
		bctrlpressed = true;
	}

	if ( _bRotateView ) {
		Vec3 spoint = get_pos_sp(x, y);
		Vec3 n = Cross(_spoint_prev, spoint);
		if ( Norm(n) > 1e-6 ) {
			if ( !bctrlpressed ) { n *= 3.0; }
			Vec3 w = ~(_mvmatrix_prev.GetRotation()) * n;
			_mvmatrix.SetRotation(_mvmatrix_prev.GetRotation() * Exp(w));
		} else {
			return;
		}
	}
	if (_bTranslateView) { float den = 30; if ( bctrlpressed ) { den = 300; } _xcam += (float)(x-_downX)/den*_mscale; ycam += (float)(_downY-y)/den*_mscale; } // translate
	if (_bZoomInOut) { float den = 30; if ( bctrlpressed ) { den = 300; } _sdepth -= (float)(_downY-y)/den*_mscale;  } // zoom in/out
	_downX = x; _downY = y;

	glutPostRedisplay();
}
 
void InitGL() 
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(640, 480);
	glutCreateWindow("trajoptim");
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glClearColor(_backgroundcolor[0],_backgroundcolor[1],_backgroundcolor[2],_backgroundcolor[3]);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHT0);
	glEnable(GL_CULL_FACE);
	glLightfv(GL_LIGHT0, GL_POSITION,_light0Position);
	glutReshapeFunc(ReshapeCallback);
	glutDisplayFunc(DisplayCallback);
	glutKeyboardFunc(KeyboardCallback);
	glutSpecialFunc(SpecialKeyboardCallback);
	glutMouseFunc(MouseCallback);
	glutMotionFunc(MotionCallback); 
	glutTimerFunc(_replay_render_freq, Timer, 0);
}
 
void LoadMySystem(int type)
{
	_pmyworld = new MyWorld;
	if ( !_pmyworld->create() ) {
		std::cout << "Error:: Falied in loading system." << std::endl;
		exit(0);
	}
	_current_time = 0.0;
	_current_frame_idx = 0;
	_mscale = 1.0;

	// set camera
	_mvmatrix = SE3(SO3(RMatrix("1 0 0; 0 0 1; 0 -1 0").GetPtr()), Vec3(0,0,-2)); // side view (eye gaze = +y)
}
 
int main(int argc, char **argv) 
{
	LoadMySystem(0);
	glutInit(&argc, argv);
	InitGL();
	glutMainLoop(); 
	return 0; 
}
 