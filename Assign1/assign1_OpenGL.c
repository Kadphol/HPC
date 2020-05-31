#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>


#ifdef __APPLE__
#include <OpenGL/gl3.h>
#define __gl_h_
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#define max(X,Y) ((X) > (Y) ? (X) : (Y))
#define screenX 1366
#define screenY 768

float r[screenX][screenY], g[screenX][screenY], b[screenX][screenY];
float matrixInput[screenX*screenY], matrixOutput[screenX*screenY];
int iter = 0,loop = 0,row,col;
int setpos,setrow,setcol;
char title[100];

void idle()
{  
    for(int i=0;i<screenX;i++)
    {
        for(int j=0;j<screenY;j++)
        {
            float ratio = 2 * (matrixInput[(i*col)+j])/255;
            r[i][j] = (float)(max(0,ratio- 1));
            b[i][j] = (float)(max(0,1 - ratio));
            g[i][j] = (float)(1 - r[i][j] - b[i][j]); 
        }
    }   
    //usleep(100000); //sleep 0.1 second
    glutPostRedisplay();         
}

void display(void)
{
    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    gluOrtho2D(0.0, screenX, 0.0, screenY);
    for(int i = 0; i < screenX; i++)
    {
        for(int j = 0; j < screenY; j++)
        {
            glColor3f(r[i][j],g[i][j],b[i][j]); 
            glBegin(GL_POINTS);
            glVertex2i (i,j);
            glEnd();
        }
    } 
    if(loop < iter)
    {
        
        for(int i = col+1; i< setpos; i++)
        {
            matrixOutput[i] = (matrixInput[(i-(col+1))] + matrixInput[(i-(col))] + matrixInput[(i-(col-1))]+
                                 matrixInput[i-1] + matrixInput[i] + matrixInput[i+1]+
                                 matrixInput[(i+(col-1))] + matrixInput[(i+(col))] + matrixInput[(i+(col+1))])/9;
        }
        for(int i = 1; i < setrow; i++)
        {
            for(int j = 1; j < setcol; j++)
            {
                matrixInput[(i*col) +j] = matrixOutput[(i*col) + j];
            }
        }
        loop++;
        if(loop < iter)
            sprintf(title,"(Time: %d) Heat Transfer Simulation",loop);
        else
            sprintf(title,"(Time: completed) Heat Transfer Simulation");

        glutSetWindowTitle(title);

    }      
    //glutSwapBuffers();
    glFlush();
}



int main(int argc, char* argv[])
{
    int rank, size;
    int i = 0, j, k;
    FILE *inputfile,*outputfile; 

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Request send_req[10],recv_req[20];

    iter = atoi(argv[3]);
    inputfile = fopen(argv[1],"r");
    fscanf(inputfile,"%d %d",&row,&col);
    while(fscanf(inputfile,"%f",&matrixInput[i++])==1);
    fclose(inputfile);
    setpos = ((row-1)* col)-1;
    setrow = row-1;
    setcol = col-1;

    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(screenX, screenY);
    glutCreateWindow("Heat Transfer Simulation");
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);

    glutDisplayFunc(display);
	glutIdleFunc(idle);
    glutMainLoop();
    
    MPI_Finalize();
    return 0;
}
