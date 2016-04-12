
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include"translation.h"

// calculates the a new mesh translated by length length
void translation2D(pMesh mesh, double lengthX, double lengthY)
{
	int i;
	pPoint ppt;
	double x,y;
	int refNew;
	for(i=0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		x = ppt->c[0] + lengthX;
		y = ppt->c[1] + lengthY;
		refNew = ppt->c[2];
	
		ppt->c[0] = x;
		ppt->c[1] = y;	
		ppt->c[2] = refNew ;	
	}
	
}

// calculates the a new mesh translated by length length
void translation3D(pMesh mesh, double lengthX, double lengthY, double lengthZ)
{
	int i;
	pPoint ppt;
	double x,y,z;
	int refNew;
	for(i=1; i <= mesh->np; i++)
	{
	
		mesh->point[i].c[0] = mesh->point[i].c[0] + lengthX;
		mesh->point[i].c[1]  = mesh->point[i].c[1] + lengthY;
		mesh->point[i].c[2] = mesh->point[i].c[2] + lengthZ;
		
	
	
	}
	
}

