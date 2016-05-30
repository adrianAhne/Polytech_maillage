
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "curvature.h"



/*
	Calculates the curvature for each point by adding up
	all the angles of the triangles around the point
	and substract this sum from 2*PI
	It creates the sol file in the current working directory	
	3 dimensional
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
					pos ++ ;
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
			
			/* Here we have the angle for each triangle around the point */
			angle +=  acos(costr) ;
						
		}
			
		mesh->sol[i] = 2.0*PI - angle ;
		free(triangle) ;
	}
	return (1) ;
	
}


/*
	Calculates the curvature for each point by adding up
	all the angles of the triangles around the point
	and substract this sum from 2*PI
	It creates the sol file in the current working directory	
	2 dimensional
*/
int courbure2D( pMesh mesh )
{
	// check tells us how many edges are connected to a vertex i
	int check = 0 , i , j ;
	int *voisin  ;
	double u1=0,u2=0,v1=0,v2=0,normu=0,normv=0,costh=0,dist1,dist2,dist;
	/* first we have to find the vertices */
	
	/* Go through all the vertices */
	for(i=1;i<= mesh->np;i++)
	{
		check = 0 ;
		
		/* MEMORY ALLOCATION */
		for(j=1;j<mesh->na;j++) 
		{
			/* here we see in which edge is the vertice i */
			if( mesh->edge[j].v[0] == i  ||  mesh->edge[j].v[1] == i )
				check ++ ;
		}
		
		voisin = (int*) calloc (  check , sizeof(int)) ;
		check =0;
		
		/* going though all the edge to find the vertices which are important */
		for(j=1;j<= mesh->na;j++)
		{
			/* here we see in which edge is the vertice i */
			if( mesh->edge[j].v[0] == i  )
			{
				/* here we stock the following point */
				voisin [1] = mesh->edge[j].v[1] ;
				check++;
			}
			
			if ( mesh->edge[j].v[1] == i )
			{
				/* here we stock the last point */
				voisin [0] = mesh->edge[j].v[0] ;
				//fprintf(stdout," pt précedent = %d \n",mesh->edge[j].v[0]);
				check++;
			}
		
		}
		
		
		if( check == 2 )
		{
			/* calculate curvature at the point */
			/* u is the vector for (pi-1,pi) and v for (pi,pi+1) */
			u1 =  mesh->point[voisin[1]].c[0] - mesh->point[i].c[0]  ;
			u2 =  mesh->point[voisin[1]].c[1] - mesh->point[i].c[1]  ;
			v1 =  mesh->point[i].c[0] - mesh->point[voisin[0]].c[0]  ;
			v2 =  mesh->point[i].c[1] - mesh->point[voisin[0]].c[1]  ;
			
			/* NORM */
			normu = sqrt(  pow(u1,2) + pow(u2,2) ) ;
			normv = sqrt(  pow(v1,2) + pow(v2,2) ) ; 
			dist1 = mesh->point[voisin[1]].c[0] - mesh->point[voisin[0]].c[0] ;
			dist2 = mesh->point[voisin[1]].c[1] - mesh->point[voisin[0]].c[1]  ;
			dist = sqrt ( pow ( dist1,2) + pow ( dist2,2));
			
			/* SCALAR PRODUCT */
			costh = ( (u1*v1) + (u2*v2) ) / ( normu * normv ) ;
			
			/* THETA */
			mesh->sol[i] = 2*sin(acos(costh));
			
			
		}
		
		// Check if 0,1 edges or more than 2 edges are connected to a vertex. If this is the case set value to -1 to receive a visual difference in the image later
		if ( check == 1 || check == 0 )
		{
			mesh->sol[i] = -1 ;			
		}
		if ( check > 2  )
		{
			mesh -> sol[i] = -1 ;
		}
	}
	
	free (voisin);
	
	return 1;

}



