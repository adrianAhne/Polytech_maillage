#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "groupfunctions.h"



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

// I don't know if it works
void center2D(Mesh *mesh, double *xc, double *yc)
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
void translation2D(Mesh *mesh, double lengthX, double lengthY)
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
void translation3D(Mesh *mesh, double lengthX, double lengthY, double lengthZ)
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
			fprintf(stdout," costr = %f \n" , acos(costr)) ;
			
			/* Here we have the angle for each triangle around the point */
			angle +=  acos(costr) ;
			
		}
		
		
	
	
		fprintf(stdout," angle = %f \n" , angle ) ;
		/* Here we save the solution for each vertice in the struct */
		
		/*if ( 2*PI - angle < 0.01 )
			 mesh->sol[i] = 0 ;*/
		//else
			mesh->sol[i] = 2.0*PI - angle ;
		
		//fprintf(stdout," solu = %f \n" , mesh->sol[i] ) ;
	}
	return (1) ;
	
}

/*** FUNCTION COURBURE ***/ 
/* In this function we will calculate the curvature for each point in a 2D Mesh.
		In order to do that we will use the gaussian curvature: this is the angular between two edge.
		PARAMETERS : the 2D mesh ( pMesh )
		RETURN : VOID 
*/
int courbure2D( pMesh mesh )
{

	int check = 0 , i , j ;
	int *voisin = 0 ;
	double u1=0,u2=0,v1=0,v2=0,normu=0,normv=0,costh=0;
	/* first we have to find the vertices */
	

	/* Go through all the vertices */
	for(i=1;i<= mesh->np;i++)
	{
		fprintf(stdout, "i = %d \n" , i ) ;
		check = 0 ;
		/* MEMORY ALLOCATION */
		for(j=1;j<mesh->na;j++)
			/* here we see in which edge is the vertice i */
			if( mesh->edge[j].v[0] == i  ||  mesh->edge[j].v[1] == i )
				check ++ ;
		voisin = (int*) malloc ( check * sizeof(int)) ;
		check =0;
		/* going though all the edge to find the vertices which are important */
		for(j=1;j<= mesh->na;j++)
		{
			/* here we see in which edge is the vertice i */
			if( mesh->edge[j].v[0] == i  )
			{
				/* here we stock the following point */
				voisin [1] = mesh->edge[j].v[1] ;
				fprintf(stdout," pt suivant = %d \n",mesh->edge[j].v[1]);
				check++;
			}
			if ( mesh->edge[j].v[1] == i )
			{
				/* here we stock the last point */
				voisin [0] = mesh->edge[j].v[0] ;
				fprintf(stdout," pt précedent = %d \n",mesh->edge[j].v[0]);
				check++;
			}
		
		}
		
		if( check == 2 )
		{
			/* calcule de la courbure au point */
			/* u is the vector for (pi-1,pi) and v for (pi,pi+1) */
			u1 =  mesh->point[voisin[0]].c[0] - mesh->point[i].c[0]  ;
			u2 =  mesh->point[voisin[0]].c[1] - mesh->point[i].c[1]  ;
			v1 =  mesh->point[i].c[0] - mesh->point[voisin[0]].c[0]  ;
			v2 =  mesh->point[i].c[1] - mesh->point[voisin[0]].c[1]  ;
			/* NORM */
			normu = sqrt(  pow(u1,2) + pow(u2,2) ) ;
			normv = sqrt(  pow(v1,2) + pow(v2,2) ) ; 
			/* SCALAR PRODUCT */
			costh = ( (u1*v1) + (u2*v2) ) / ( normu * normv ) ;
			/* THETA */
			mesh->sol[i] = acos(costh) ;
			
			
			
			/* pondération par les longueurs (li-1,li) et (li,li+1)*/
		}
		if ( check == 1 )
		{
			fprintf(stdout, " Point pas à l'intersection de 2 edges mais fait partie de l'un d'entre eux \n" );
			mesh->sol[i] = 0.0 ;
			fprintf(stdout," check = %d \n" ,check );
			
		}
		if ( check > 2 )
		{
			fprintf(stdout ," PROBLEME point à l'intersection de plusieurs edges \n" );	
			mesh -> sol[i] = 0.0 ;
			
		}
		free (voisin);
	}
	return 1;

}

void normalesOfTriangles(Mesh *mesh)
{

	mesh->triaNorm = (pTriaNorm)calloc(mesh->nt+1, sizeof(TriaNorm));
	assert(mesh->triaNorm);

	mesh->Normal = (pNormal)calloc(mesh->np+1, sizeof(Normal));
	assert(mesh->Normal);

	int i, j, k;
	pTria currentTria;
	Point P1,P2,P3,N;
	double weight, norm;

	// Loop over all triangles calculating the normales, the area and store data into mesh->triaNorm
	for(i=1; i <= mesh->nt; i++)
	{
		currentTria = &mesh->tria[i];
		
		/*
		In this part, I'm calculating the normal of the face of the triangle
		After, we'll average each normal
		For this technique, we'll have to calculate the weighting of each face (= surface ratio)

		Abstract : normal of 1 triangle (composed of 3 points P1, P2 and P3) :

		The cross product of two sides of the triangle equals the surface normal. 
		So, if V = P2 - P1 and W = P3 - P1, and N is the surface normal, then:

		Nx=(Vy∗Wz)−(Vz∗Wy)=((P2y-P1y)*(P3z-P1z)-(P2z-P1z)*(P3y-P1y))
		Ny=(Vz∗Wx)−(Vx∗Wz)=((P2z-P1z)*(P3x-P1x)-(P2x-P1x)*(P3z-P1z))
		Nz=(Vx∗Wy)−(Vy∗Wx)=((P2x-P1x)*(P3y-P1y)-(P2y-P1y)*(P3x-P1x))

		Weight=AreaOfTheTriangle=0.5*sqrt((Vy*Wz-Vz*Wy)^2+(Vz*Wx-Vx*Wz)^2+(Vx*Wy-Vx*Wx)^2)

		Translation in code :
		*/

		P1 = mesh->point[(currentTria->v[0])];
		P2 = mesh->point[(currentTria->v[1])];
		P3 = mesh->point[(currentTria->v[2])];
		
		// Normal of the triangle
		N.c[0] = (P2.c[1]-P1.c[1])*(P3.c[2]-P1.c[2]) - (P2.c[2]-P1.c[2])*(P3.c[1]-P1.c[1]);
		N.c[1] = (P2.c[2]-P1.c[2])*(P3.c[0]-P1.c[0]) - (P2.c[0]-P1.c[0])*(P3.c[2]-P1.c[2]);
		N.c[2] = (P2.c[0]-P1.c[0])*(P3.c[1]-P1.c[1]) - (P2.c[1]-P1.c[1])*(P3.c[0]-P1.c[0]);

		// Weight of the triangle (let's say it the area)
		weight = 0.5 * sqrt( pow(N.c[0],2)+pow(N.c[1],2)+pow(N.c[2],2) );
		
		// Normale is calculated at point P1 for the triangle; save both points defining the normale
		(mesh->triaNorm[i]).n[0] = N.c[0];
		(mesh->triaNorm[i]).n[1] = N.c[1];
		(mesh->triaNorm[i]).n[2] = N.c[2];
		printf("%f %f %f\n", N.c[0],N.c[1],N.c[2]);
		(mesh->triaNorm[i]).weight = weight;
		

	}


	// REMARQUES
	// Faire plutot une boucle sur les triangles
	// Calculer l'aire à chaque boucle après
	//
	// A la fin du calcul de normales, normaliser chaque normale pour chaque point par sa norme

	
	// Loop over all points and calculate weighted normales for each point by using neighbour triangles
	// WARNING : Not optimized and not normalized normales ! 
	/*
	for (i = 1; i <= mesh->np; i++)
	{
		mesh->point[i].n[0] = 0;
		mesh->point[i].n[1] = 0;
		mesh->point[i].n[2] = 0;

		for (j = 0; j <= mesh->nt; j++)
		{
			currentTria = &mesh->tria[j];
			for (k = 0; k < 3; k++)
			{
				if (currentTria->v[k] == i)
				{

					mesh->point[i].n[0] += mesh->triaNorm[j].weight * (mesh->triaNorm[j].n[0]);
					mesh->point[i].n[1] += mesh->triaNorm[j].weight * (mesh->triaNorm[j].n[1]);
					mesh->point[i].n[2] += mesh->triaNorm[j].weight * (mesh->triaNorm[j].n[2]);
				}
			}
		}
		mesh->point[i].n[0] = mesh->point[i].n[0] / sqrt( pow(mesh->point[i].n[0],2) + pow(mesh->point[i].n[1],2) + pow(mesh->point[i].n[2],2) );
		mesh->point[i].n[1] = mesh->point[i].n[1] / sqrt( pow(mesh->point[i].n[0],2) + pow(mesh->point[i].n[1],2) + pow(mesh->point[i].n[2],2) );
		mesh->point[i].n[2] = mesh->point[i].n[2] / sqrt( pow(mesh->point[i].n[0],2) + pow(mesh->point[i].n[1],2) + pow(mesh->point[i].n[2],2) );

	}
	*/

	for (i = 1; i < mesh->nt ; i++) {
		for (j = 0; j < 3 ; j++) {
			for (k = 0; k < 3; k++)
			{
				mesh->Normal[mesh->tria[i].v[j]].n[k] += mesh->triaNorm[i].n[k];
			}
						
		}
	}

	for (i = 1; i < mesh->np; ++i)
	{
		norm = sqrt( pow(mesh->Normal[i].n[0],2) + pow(mesh->Normal[i].n[1],2) + pow(mesh->Normal[i].n[2],2) );
		for (j=0;j<3;j++) {
			mesh->Normal[i].n[j] *= 1 / norm;
		}
	}

	mesh->nn = mesh->np;

	

}



