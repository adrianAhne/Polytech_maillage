#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"



void rotation2D(Mesh *mesh, double angle)
{
	pPoint ppt;
	int i;

	//printf("Centres : %f %f\n", xc, yc);
	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		double x,y;
		// Check if 2D or 3D
		
		x = ppt->c[0] * cos(angle) - ppt->c[1] * sin(angle);
		y = ppt->c[1] * cos(angle) + ppt->c[0] * sin(angle);

		ppt->c[0] = x;
		ppt->c[1] = y;
	}
}


void rotation3D(Mesh *mesh, double angleX, double angleY, double angleZ)
{
	pPoint ppt;
	int i;

	double x,y,z;
	
	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		
		
		// Si l'angle est sufisamment grand pour que la rotation puisse se faire
		if (fabsf(angleX) > 1e-10)
		{
			y = ppt->c[1] * cos(angleX) - ppt->c[2] * sin(angleX);
			z = ppt->c[1] * sin(angleX) + ppt->c[2] * cos(angleX);
			ppt->c[1] = y;
			ppt->c[2] = z;
		}

		if (fabsf(angleY) > 1e-10)
		{
			x = ppt->c[2] * sin(angleY) + ppt->c[0] * cos(angleY);
			z = ppt->c[2] * cos(angleY) - ppt->c[0] * sin(angleY);
			ppt->c[0] = x;
			ppt->c[2] = z;
		}

		if (fabsf(angleZ) > 1e-10)
		{
			x = ppt->c[0] * cos(angleZ) - ppt->c[1] * sin(angleZ);
			y = ppt->c[0] * sin(angleZ) + ppt->c[1] * cos(angleZ);
			ppt->c[0] = x;
			ppt->c[1] = y;
		}
		
	}
}

