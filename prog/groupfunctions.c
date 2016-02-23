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

	//printf("Centres : %f %f\n", xc, yc);
	for (i = 0; i <= mesh->np; i++)
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

<<<<<<< HEAD

=======
>>>>>>> 82c2e7b563274d93f549b1218f7e8c6ed39f5b39
void rotation3D(Mesh *mesh, float angleX, float angleY, float angleZ)
{
	pPoint ppt;
	int i;

	float x,y,z;
	
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


/*** Change 2D to 3D ***/
/* Parameter : a mesh with dim == 2 
	This function will add a square to the tab of points for the third dimension
*/
void Change2Dto3D( pMesh Mesh )
{
	int i;
	for(i=0;i<=Mesh -> np;i++)
	{
		Mesh->point[i].c[3] = Mesh->point[i].c[2] ;
		Mesh->point[i].c[2] = 0.0 ;
	}
}

<<<<<<< HEAD
=======

>>>>>>> 081dae967ff65a11838a62a4187303f0399c1a75
/*** FUNCTION SUPERPOSITION ***/
/*	Parameters : 3 meshs: the 2 meshs to combine + the combinaison of the 2
	return : 1 if it works, 0 if not
	This function will join a mesh to another 
	Before this function The function loadMesh has to be done on the two meshs to combine, 
	in order to fill the number of Vertices, Triangles and edges of each one.
	After that we will go through the file of Mesh1 as a text file and find where is the last point
*/

int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) 
{
	
	/* Variables */
	int i,j;
	
	/* First we put the number of vertices, triangles and edges of the final mesh */ 
	Mesh_final->nt = Mesh1->nt + Mesh2->nt 	;
	Mesh_final->np = Mesh1->np + Mesh2->np 	;
	Mesh_final->na = Mesh1->na + Mesh2->na 	;
	Mesh_final->nr = Mesh1->nr + Mesh2->nr 	;

	/* Fill the dimension ( we take the more important dimension) and the mark */

	if ( Mesh1->dim == 3 && Mesh2->dim == 2 )
	{
		Mesh_final-> dim = 3 ;
		Change2Dto3D(Mesh2);
	}
	if (Mesh2->dim == 3 && Mesh1->dim == 2)
	{
		Mesh_final-> dim = 3 ;
		Change2Dto3D(Mesh1) ;
	}
	if (Mesh1->dim == 3 && Mesh2->dim == 3)
		Mesh_final-> dim = 3 ;
	
	if (Mesh1->dim == 2 && Mesh2->dim == 2)
		Mesh_final-> dim = 2 ;
	
	Mesh_final->mark = Mesh1->mark ;
	
	
	/* Version of the mesh wanted */
	Mesh_final->ver = 1 ;
	

	/* Memory allocation */

	Mesh_final->point = (pPoint)calloc(Mesh_final->np,sizeof(Point));
  assert(Mesh_final->point);
	Mesh_final->sol = (double*)calloc(3*Mesh_final->np, sizeof(double));
	assert(Mesh_final->sol);
  if ( Mesh_final->na ){
    Mesh_final->edge = (pEdge)calloc(Mesh_final->na,sizeof(Edge));
    assert(Mesh_final->edge);
  }
  if ( Mesh_final->nt ) {
    Mesh_final->tria = (pTria)calloc(Mesh_final->nt,sizeof(Tria));
    assert(Mesh_final->tria);
  }
  
  
	
	/* Now we fill the tab of points, vertices and edges */

	/* For the tab of points we first input in the new tab the vertices of mesh1. 
		 After that we eliminate the last char of the tab ("\0") and input the vertices of the mesh2
	*/  
	
	
		int tab = 0 ;
	
		for(i=0;i<=Mesh1->np;i++)
		{
			Mesh_final->point[tab] = Mesh1->point[i];
			tab ++ ;
		}
		
	/*points of mesh2*/
	
		for(j=0;j<= Mesh2->np;j++)
		{
			Mesh_final->point[tab] = Mesh2->point[j+1];
			tab ++ ;
		}
	
	tab = 0 ;
	/* Triangles of mesh1 */
	
		for(i=0;i<=Mesh1->nt;i++)
		{
			Mesh_final->tria[tab] = Mesh1->tria[i];
			tab++ ;
		}

	/*Triangles of mesh2*/
	
		for(j=0;j<= Mesh2->nt;j++)
		{
			Mesh_final->tria[tab].v[0] = (Mesh2->tria[j+1].v[0] + Mesh1->np );
			Mesh_final->tria[tab].v[1] = (Mesh2->tria[j+1].v[1] + Mesh1->np );
			Mesh_final->tria[tab].v[2] = (Mesh2->tria[j+1].v[2] + Mesh1->np );
			tab ++ ;
		}
	return(1);
}

// calculates the a new mesh translated by length length
void translation2D(Mesh *mesh, float lengthX, float lengthY)
{
	int i;
	pPoint ppt;
	float x,y;
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
void translation3D(Mesh *mesh, float lengthX, float lengthY, float lengthZ)
{
	int i;
	pPoint ppt;
	float x,y,z;
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
