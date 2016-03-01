#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "groupfunctions.h"


/*** FUNCTION ROTATION ***/
void rotation2D(pMesh mesh, double angle)
{
	pPoint ppt;
	int i;

	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		ppt->c[0] = ppt->c[0] * cos(angle) - ppt->c[1] * sin(angle);
		ppt->c[1] = ppt->c[1] * cos(angle) + ppt->c[0] * sin(angle);
    
	}
}

void rotation3D(pMesh mesh, double angleX, double angleY, double angleZ)
{
	pPoint ppt;
	int i;
	
	for (i = 0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
		
		// Si l'angle est sufisamment grand pour que la rotation puisse se faire
		if (fabs(angleX) > 1e-10)
		{
			ppt->c[1] = ppt->c[1] * cos(angleX) - ppt->c[2] * sin(angleX);
			ppt->c[2] = ppt->c[1] * sin(angleX) + ppt->c[2] * cos(angleX);
		
		}

		if (fabs(angleY) > 1e-10)
		{
			ppt->c[0] = ppt->c[2] * sin(angleY) + ppt->c[0] * cos(angleY);
			ppt->c[2] = ppt->c[2] * cos(angleY) - ppt->c[0] * sin(angleY);

		}

		if (fabs(angleZ) > 1e-10)
		{
			ppt->c[0] = ppt->c[0] * cos(angleZ) - ppt->c[1] * sin(angleZ);
			ppt->c[1] = ppt->c[0] * sin(angleZ) + ppt->c[1] * cos(angleZ);
		}
		
	}
}

/* The mesh center c is the centre of the bounding box [xmin,xmax]*[ymin,ymax]*[zmin,zmax] */
int center(pMesh mesh, double *c)//ici j'ai mis c (au lieu que xc,yc) pour traiter les cas 2d et 3d dans la même fonction
{
	pPoint  ppt;
	int     i,l;
  double   min[3],max[3];
  
  /*compute bounding box*/
  for (i=0; i<mesh->dim; i++)
  {
    min[i] =   FLT_MAX;
    max[i] =  -FLT_MAX;
    
  }
	for (i = 1; i <= mesh->np; ++i)
	{
    ppt = &mesh->point[i];
    for (l=0; l<mesh->dim; l++)
    {
      max[l] = D_MAX(max[l],ppt->c[l]);
      min[l] = D_MIN(min[l],ppt->c[l]);
    }
  }
  if(mesh->dim==3)
  fprintf(stdout,"  %%%% Bounding box  [%f,%f] [%f,%f] [%f,%f] \n",min[0],max[0],min[1],max[1],min[2],max[2]);
  else
   fprintf(stdout,"  %%%% Bounding box  [%f,%f] [%f,%f] \n",min[0],max[0],min[1],max[1]);
  
  /*compute center of the bounding box*/
  for (l=0; l<mesh->dim; l++) c[l] = 0.5*(min[l]+max[l]);
  
  if(mesh->dim==3)
    fprintf(stdout,"  %%%% Center  %f %f %f \n",c[0],c[1],c[2]);
  else
    fprintf(stdout,"  %%%% Center  %f %f \n",c[0],c[1]);
  
  
    return 1;
}

/******************************/
/*** FUNCTION SUPERPOSITION ***/
/******************************/
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

	Mesh_final->dim = Mesh1->dim ;
	Mesh_final->mark = Mesh1->mark ;
	
	
	/* Version of the mesh wanted */
	Mesh_final->ver = 1 ;
	

	/* Memory allocation */

	Mesh_final->point = (pPoint)calloc(Mesh_final->np,sizeof(Point));
  assert(Mesh_final->point);
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
	
	/*points of mesh1*/
		
    for(i=1;i<=Mesh1->np;i++)
		{
			Mesh_final->point[i] = Mesh1->point[i];
		}
	
	/*points of mesh2*/
  
    for(j=1;j<= Mesh2->np;j++)
		{
			Mesh_final->point[i+j] = Mesh2->point[j];
		}
	
	/* Triangles of mesh1 */
  
    for(i=1;i<=Mesh1->nt;i++)
		{
			Mesh_final->tria[i] = Mesh1->tria[i];
		}

	/*Triangles of mesh2*/
	
		for(j=1;j<= Mesh2->nt;j++)
		{
			Mesh_final->tria[i+j].v[0] = (Mesh2->tria[j+1].v[0])+((Mesh1->np) + 1);
			Mesh_final->tria[i+j].v[1] = (Mesh2->tria[j+1].v[1])+((Mesh1->np) + 1);
			Mesh_final->tria[i+j].v[2] = (Mesh2->tria[j+1].v[2])+((Mesh1->np) + 1);
		}
	return(1);
}

/****************************/
/*** FUNCTION TRANSLATION ***/
/****************************/

// calculates a new mesh translated by length length
void translation2D(pMesh mesh, double lengthX, double lengthY)
{
	int i;
	pPoint ppt;
  
	for(i=0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
	  ppt->c[0] += lengthX;
		ppt->c[1] += lengthY;
	}
}

// calculates a new mesh translated by length length
void translation3D(Mesh *mesh, double lengthX, double lengthY, double lengthZ)
{
	int i;
	pPoint ppt;

	for(i=0; i <= mesh->np; i++)
	{
		ppt = &mesh->point[i];
	
    ppt->c[0] += lengthX;
		ppt->c[1] += lengthY;
		ppt->c[2] += lengthZ;
	}
	
}


/*** FUNCTION COURBURE ***/
/* 
	This functions will calculate the curve for each point 
	For that we will add up all the angle of the triangles around the point and make the difference to 2PI 
	Parameters : a mesh 
	return : it create the sol file in the current repository
						it return if it works ,  if not 
*/
int courbure3D(pMesh mesh )
{
	
	/* DEFINITION */
	double angle,distAB,distAC,costr ;
	int i,j,point1,point2, THEpoint ;
	pTria	triangle ;
	int pos = 0 ;
	double AB[3],AC[3] ;
	

	
	
	/* CALCUL */
	for (i=1 ; i<= mesh->np ; i++ )
	{
	
		angle = 0.000000000000 ;

		
		/* for each point of the mesh we stock in a tab the triangles with the point */ 
		for (j=1 ; j <= mesh -> nt ; j++)
		{		
			
			if ( mesh -> tria[j].v[0] == i || mesh -> tria[j].v[1] == i || mesh -> tria[j].v[2] == i )
			{
					/*triangle[pos] = mesh->tria[j] ;*/
					pos ++ ;
			}
		}
			/* MEMORY ALLOCATION */
	
		triangle = ( pTria ) calloc ( (pos+1) , sizeof(Tria) ) ;
		pos = 0 ;
		
		for (j=1 ; j <= mesh -> nt ; j++)
		{		
			
			if ( mesh -> tria[j].v[0] == i || mesh -> tria[j].v[1] == i || mesh -> tria[j].v[2] == i )
			{		
					
					triangle[pos] = mesh->tria[j] ;
					pos ++ ;
			}
			
		}
		//fprintf(stdout," pos = %d \n" , pos );
		
		/* first we calcule the cosinus of each angle around the point */
		for (j=0 ; j < pos ; j++ )
		{
			/* we determine which point is the current point */
			if  (   triangle[j].v[0] == i )
			{
				THEpoint = 0;
				point1 = 1 ;
				point2 = 2 ;
			}
			if  (   triangle[j].v[1] == i )
			{
				THEpoint = 1;
				point1 = 0 ;
				point2 = 2 ;
			}
			if  (   triangle[j].v[2] == i )
			{
				THEpoint = 2;	
				point1 = 0 ;
				point2 = 1 ;
			}
		
			/* now we calculate the scalar product */
			
			
			/*
			fprintf(stdout," coordonnées de A = ( %f , %f , %f ) \n ", mesh->point[i].c[0] , mesh->point[i].c[1] , mesh->point[i].c[2] ) ;
			fprintf(stdout," coordonnées de B %d = ( %f , %f , %f ) \n ", triangle[j].v[point1] , mesh->point[triangle[j].v[point1]].c[0] , mesh->point[triangle[j].v[point1]].c[1] , mesh->point[triangle[j].v[point1]].c[2] ) ;
			fprintf(stdout," coordonnées de C %d = ( %f , %f , %f ) \n ", triangle[j].v[point2], mesh->point[triangle[j].v[point2]].c[0] , mesh->point[triangle[j].v[point2]].c[1] , mesh->point[triangle[j].v[point2]].c[2] ) ;
		*/
			AB[0] =  mesh->point[triangle[j].v[point1]].c[0] - mesh->point[i].c[0]  ;
			AB[1] =  mesh->point[triangle[j].v[point1]].c[1] - mesh->point[i].c[1]  ;
			AB[2] =  mesh->point[triangle[j].v[point1]].c[2] - mesh->point[i].c[2]  ;
			
			AC[0] = mesh->point[triangle[j].v[point2]].c[0] - mesh->point[i].c[0]  ;
			AC[1] = mesh->point[triangle[j].v[point2]].c[1] - mesh->point[i].c[1]  ;
			AC[2] = mesh->point[triangle[j].v[point2]].c[2] - mesh->point[i].c[2]  ;
			
			distAB = sqrt(  pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2)  ) ;
			distAC = sqrt(  pow(AC[0],2) + pow(AC[1],2) + pow(AC[2],2)  ) ; 
			
			
			if ( distAB < 0 || distAC < 0 )
				fprintf(stdout, "PROBLEME ! \n");
			costr = ( AB[0] * AC [0] + AB[1] * AC [1] + AB[2] * AC[2] ) / (  distAB * distAC ) ;
			//fprintf(stdout," costr = %f \n" , costr) ;
			
			/* Here we have the angle for each triangle around the point */
			angle +=  acos(costr) ;
			
		}
		
		
	
	
		//fprintf(stdout," angle = %f \n" , angle ) ;
		/* Here we save the solution for each vertice in the struct */
		
		/*if ( 2*PI - angle < 0.01 )
			 mesh->sol[i] = 0 ;*/
		//else
			mesh->sol[i] = 2.0*PI - angle ;
		
		//fprintf(stdout," solu = %f \n" , mesh->sol[i] ) ;
	}
	return (1) ;
	
}

void normalesOfTriangles(Mesh *mesh)
{
	mesh->triaNorm = (pTriaNorm)calloc(mesh->nt+1, sizeof(TriaNorm));
	assert(mesh->triaNorm);
	
	int i, j;
	pTria currentTria;
	for(i=0; i <= mesh->nt; i++)
	{
		currentTria = &mesh->tria[i];
		for(j=0; j < 3; j++)
		{
			(mesh->triaNorm[i]).point[j] = &mesh->point[(currentTria->v[j])];
		}
		// the line is faulse, check formula
		(mesh->triaNorm[i]).norm[0]->v[0] = ( mesh->point[(currentTria->v[0])].c[0] - mesh->point[(currentTria->v[1])].c[0] ) *( mesh->point[(currentTria->v[0])].c[0] - mesh->point[(currentTria->v[2])].c[0] );
			
		
		

	}


	

}









	









	


