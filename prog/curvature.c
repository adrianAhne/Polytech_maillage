
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
int curvature3D(pMesh mesh)
{
	/* For each triangle we calculate the 3 angles for each 3 points of the triangle 
		we fill a tab of point with the sum of this angles */
	int i,j ;
	double AB[3] , AC[3] , distAB , distAC ;


	/* First we declare the tab and init it at zero */
	double* tab = (double*) calloc ( (mesh->np)+1 , sizeof(double)) ;

	/* For the triangle */
	for ( i=1 ; i<= mesh->nt ; i++ )
	{
		/* go through the points of the triangle */
		for(j=0 ; j<=2 ; j++)
		{
			/* calculation of the angle */
			if ( j == 0)
			{
				AB[0] =	mesh->point[mesh->tria[i].v[1]].c[0] - mesh->point[mesh->tria[i].v[0]].c[0] ;
				AB[1] = mesh->point[mesh->tria[i].v[1]].c[1] - mesh->point[mesh->tria[i].v[0]].c[1] ;
				AB[2] = mesh->point[mesh->tria[i].v[1]].c[2] - mesh->point[mesh->tria[i].v[0]].c[2] ;

				AC[0] = mesh->point[mesh->tria[i].v[2]].c[0] - mesh->point[mesh->tria[i].v[0]].c[0] ;
				AC[1] = mesh->point[mesh->tria[i].v[2]].c[1] - mesh->point[mesh->tria[i].v[0]].c[1] ;
				AC[2] = mesh->point[mesh->tria[i].v[2]].c[2] - mesh->point[mesh->tria[i].v[0]].c[2] ;

				distAB = sqrt(  pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2)  ) ;
				distAC = sqrt(  pow(AC[0],2) + pow(AC[1],2) + pow(AC[2],2)  ) ; 

				if ( distAB < 0 || distAC < 0 )
					fprintf(stdout, "PROBLEME ! \n");
					
				tab[mesh->tria[i].v[0]] += acos(( AB[0] * AC [0] + AB[1] * AC [1] + AB[2] * AC[2] ) / (  distAB * distAC ) );
				
				/* Here we have the angle for each triangle around the point */
				
			}
			if ( j == 1 )
			{
				AB[0] =	mesh->point[mesh->tria[i].v[0]].c[0] - mesh->point[mesh->tria[i].v[1]].c[0] ;
				AB[1] = mesh->point[mesh->tria[i].v[0]].c[1] - mesh->point[mesh->tria[i].v[1]].c[1] ;
				AB[2] = mesh->point[mesh->tria[i].v[0]].c[2] - mesh->point[mesh->tria[i].v[1]].c[2] ;

				AC[0] = mesh->point[mesh->tria[i].v[2]].c[0] - mesh->point[mesh->tria[i].v[1]].c[0] ;
				AC[1] = mesh->point[mesh->tria[i].v[2]].c[1] - mesh->point[mesh->tria[i].v[1]].c[1] ;
				AC[2] = mesh->point[mesh->tria[i].v[2]].c[2] - mesh->point[mesh->tria[i].v[1]].c[2] ;

				distAB = sqrt(  pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2)  ) ;
				distAC = sqrt(  pow(AC[0],2) + pow(AC[1],2) + pow(AC[2],2)  ) ; 

				if ( distAB < 0 || distAC < 0 )
					fprintf(stdout, "PROBLEME ! \n");
					
				tab[mesh->tria[i].v[1]] += acos(( AB[0] * AC [0] + AB[1] * AC [1] + AB[2] * AC[2] ) / (  distAB * distAC ) );
				
				/* Here we have the angle for each triangle around the point */
				
			}
			if ( j == 2 )
			{
				AB[0] =	mesh->point[mesh->tria[i].v[1]].c[0] - mesh->point[mesh->tria[i].v[2]].c[0] ;
				AB[1] = mesh->point[mesh->tria[i].v[1]].c[1] - mesh->point[mesh->tria[i].v[2]].c[1] ;
				AB[2] = mesh->point[mesh->tria[i].v[1]].c[2] - mesh->point[mesh->tria[i].v[2]].c[2] ;

				AC[0] = mesh->point[mesh->tria[i].v[0]].c[0] - mesh->point[mesh->tria[i].v[2]].c[0] ;
				AC[1] = mesh->point[mesh->tria[i].v[0]].c[1] - mesh->point[mesh->tria[i].v[2]].c[1] ;
				AC[2] = mesh->point[mesh->tria[i].v[0]].c[2] - mesh->point[mesh->tria[i].v[2]].c[2] ;

				distAB = sqrt(  pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2)  ) ;
				distAC = sqrt(  pow(AC[0],2) + pow(AC[1],2) + pow(AC[2],2)  ) ; 

				if ( distAB < 0 || distAC < 0 )
					fprintf(stdout, "PROBLEME ! \n");
					
				tab[mesh->tria[i].v[2]] += acos(( AB[0] * AC [0] + AB[1] * AC [1] + AB[2] * AC[2] ) / (  distAB * distAC ) );
				
				/* Here we have the angle for each triangle around the point */
				
			}
		}
		
	}
	for ( i=1 ; i<=mesh->np ; i++ )
	{
		
		mesh->sol[i] = 2.0*PI - tab[i] ;
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
int curvature2D(pMesh mesh)
{
	/* Here we go through all the edges and we find the 3 points */
	int i, p_0, p_i, p_2,j, test=0 ;
	double AB[2],AC[2], normAB , normAC, BC[2] , normBC;
	double* tab = (double*) calloc ( mesh->np +1 , sizeof(double)) ;
	for(i=1 ;  i <= mesh->na ; i++ ) 
	{
		p_0 = mesh->edge[i].v[0];
		p_i = mesh->edge[i].v[1];
		j = 0 ;
		test = 0 ;
		while ( !test )
		{
			if ( j != i )
			{
				if( mesh->edge[j].v[0] == p_i )
				{
					test = 1 ;
					p_2 = mesh->edge[j].v[1] ;
				}
				if( mesh->edge[j].v[1] == p_i)
				{
					test = 1 ;
					p_2 = mesh->edge[j].v[0] ;
				}
			}
			j++ ;
		} 
		
		
		/* we have the 3 points */
		AB[0] = mesh->point[p_i].c[0] - mesh->point[p_0].c[0] ; 
		AB[1] = mesh->point[p_i].c[1] - mesh->point[p_0].c[1] ; 
		
		
		AC[0] = mesh->point[p_i].c[0] - mesh->point[p_2].c[0] ; 
		AC[1] = mesh->point[p_i].c[1] - mesh->point[p_2].c[1] ; 
		

		BC[0] = mesh->point[p_0].c[0] - mesh->point[p_2].c[0] ; 
		BC[1] = mesh->point[p_0].c[1] - mesh->point[p_2].c[1] ; 
		
		normAB = sqrt( pow(AB[0],2) + pow(AB[1],2)); 
		normAC = sqrt( pow(AC[0],2) + pow(AC[1],2) );
		normBC = sqrt( pow(BC[0],2) + pow(BC[1],2) );
		/* now we need the area of the triangle to do the calculation */		
		double weight,denom, p ;
		p = (normAB + normBC + normAC ) * 0.5 ;
		
		
		
		weight = sqrt (p*(p-normAB)*(p-normAC)*(p-normBC));
		denom = 1.0 *  normAB * normBC * normAC  ;
		mesh->sol[p_i] = 4.0*weight/denom  ;

	}
	return (1);
}

