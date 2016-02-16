#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "groupfunctions.h"



void rotation2D(Mesh *mesh, float angle)
{
	pPoint ppt;
	int i;
	float xc, yc;

	printf("Centres : %f %f\n", xc, yc);
	for (i = 0; i <= mesh->np; ++i)
	{
		ppt = &mesh->point[i];
		float x,y;
		// Check if 2D or 3D
		

		x = ppt->c[0] * cos(angle) - ppt->c[1] * sin(angle);
		y = ppt->c[1] * cos(angle) + ppt->c[0] * sin(angle);

		ppt->c[0] = x;
		ppt->c[1] = y;
	}
}

// I don't know if it works
void center2D(Mesh *mesh, float *xc, float *yc)
{
	pPoint ppt;
	int i;
	for (i = 0; i <= mesh->np; ++i)
	{
		ppt = &mesh->point[i];
		*xc += ppt->c[0];
		*yc += ppt->c[1];
	}
	*xc = *xc / mesh->np;
	*yc = *yc / mesh->np;
}

