#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include"superposition.h"

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

