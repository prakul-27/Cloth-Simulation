#include "Cloth.h"
#include <GL/glut.h>
#include <FreeImage.h>
#include <unistd.h> 


void DisplayCallback(void) {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	float center = (MESH_SIZE - 1) * STRUCTURAL_NATURAL_LENGTH / 2.0f;
	gluLookAt(center, center, center * 3, center, center, 0.0f, 0.0f, 1.0f, 0.0f); // Look at the center of the cloth
	static Cloth cloth;
	cloth.Simulate();
	cloth.Draw();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0f, 1.333f, 1.0f, center * 20);

	glutSwapBuffers();
}

int main(int argc, char *argv[]) {    
    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(1024, 768);
	glutCreateWindow("Cloth Simulation");
	glutDisplayFunc(DisplayCallback);
	glutIdleFunc(DisplayCallback);

	glEnable(GL_DEPTH_TEST);
	//glShadeModel(GL_SMOOTH);

	//SetTextures("Fabric_Texture.jpg");
	//SetLightSource();
	//SetMaterial();

	glutMainLoop();     
}