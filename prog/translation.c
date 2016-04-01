
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
	for(i=0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		x = ppt->c[0] + lengthX;
		y = ppt->c[1] + lengthY;
		z = ppt->c[2] + lengthZ;
		refNew = ppt->c[3];
	
		ppt->c[0] = x;
		ppt->c[1] = y;	
		ppt->c[2] = z;		
		ppt->c[3] = refNew ;	
	}
	
}

