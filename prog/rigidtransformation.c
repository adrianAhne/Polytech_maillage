#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "rigidtransformation.h"

/*
	File contains functions to translate or rotate a mesh and to realise the superposition between two meshes into a new one.
	
*/


/* translates every point in the mesh by lengthX in x-direction and lengthY in ydirection in the 2 dimensional case */
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


/* translates every point in the mesh by lengthX in x-direction and lengthY in y-direction and lengthZ in z-direction in the 3 dimensional case*/
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



/* rotation via the angle angle in the two dimension case */
void rotation2D(Mesh *mesh, double angle)
{
	pPoint ppt;
	int i;
	double x,y;

	// rotation over all points
	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		
		x = ppt->c[0] * cos(angle) - ppt->c[1] * sin(angle);
		y = ppt->c[1] * cos(angle) + ppt->c[0] * sin(angle);

		ppt->c[0] = x;
		ppt->c[1] = y;
	}
}

/* rotation via the angles angleX, angleY, angleZ in the three dimension case */
void rotation3D(Mesh *mesh, double angleX, double angleY, double angleZ)
{
	pPoint ppt;
	int i;

	double x,y,z;
	
	// rotation over all points
	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		
		// Check if the angle is big enough to realise the rotation
		if (fabs(angleX) > 1e-10)
		{
			y = ppt->c[1] * cos(angleX) - ppt->c[2] * sin(angleX);
			z = ppt->c[1] * sin(angleX) + ppt->c[2] * cos(angleX);
			ppt->c[1] = y;
			ppt->c[2] = z;
		}

		if (fabs(angleY) > 1e-10)
		{
			x = ppt->c[2] * sin(angleY) + ppt->c[0] * cos(angleY);
			z = ppt->c[2] * cos(angleY) - ppt->c[0] * sin(angleY);
			ppt->c[0] = x;
			ppt->c[2] = z;
		}

		if (fabs(angleZ) > 1e-10)
		{
			x = ppt->c[0] * cos(angleZ) - ppt->c[1] * sin(angleZ);
			y = ppt->c[0] * sin(angleZ) + ppt->c[1] * cos(angleZ);
			ppt->c[0] = x;
			ppt->c[1] = y;
		}
		
	}
}




/* Changes from 2D to 3D by setting the z-coordinate to 0 */
void Change2Dto3D( pMesh Mesh )
{
	int i;
	for(i=0;i<=Mesh -> np;i++)
		Mesh->point[i].c[2] = 0.0 ;
}



/* Combines all information (points, edges, triangles,...) of two meshes in a new mesh*/
int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) 
{
	
	/* Variables */
	int i,j;
	
	/* First we put the number of vertices, triangles and edges to the final mesh */ 
	Mesh_final->nt = Mesh1->nt + Mesh2->nt 	;
	Mesh_final->np = Mesh1->np + Mesh2->np 	;
	Mesh_final->na = Mesh1->na + Mesh2->na 	;
	Mesh_final->nr = Mesh1->nr + Mesh2->nr 	;

	/* Transfer the dimension ( we take the more important dimension) and the mark */
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
	
	
	/* Version of the mesh */
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
  
  
	

	
		/* Transfer first the points of Mesh1 and then of Mesh2 */
		int tab = 0 ;
	
		for(i=0; i<=Mesh1->np; i++)
		{
			Mesh_final->point[tab] = Mesh1->point[i];
			tab ++ ;
		}
		
		for(j=0;j<= Mesh2->np;j++)
		{
			Mesh_final->point[tab] = Mesh2->point[j+1];
			tab ++ ;
		}
	
	
		/* Transfer first the triangles of Mesh1 and then of Mesh2 */
		tab = 0 ;
	
		for(i=0;i<=Mesh1->nt;i++)
		{
			Mesh_final->tria[tab] = Mesh1->tria[i];
			tab++ ;
		}
	
		for(j=0;j<= Mesh2->nt;j++)
		{
			Mesh_final->tria[tab].v[0] = (Mesh2->tria[j+1].v[0] + Mesh1->np );
			Mesh_final->tria[tab].v[1] = (Mesh2->tria[j+1].v[1] + Mesh1->np );
			Mesh_final->tria[tab].v[2] = (Mesh2->tria[j+1].v[2] + Mesh1->np );
			tab ++ ;
		}
		
	return(1);
}
