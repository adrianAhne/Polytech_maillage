/*	
		CONSTRUCTION BOULE
		
		Jean-Tupac Quiroga
		Adrian Ahne
		Alexandre POULAIN
		
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include <time.h>
#include "ball.h"



/* Searches the ball/triangles around the given point and saves them in list
   Brute force method, walking trough all triangles and check if points is 
   attached or not
   triangles are saved like this: list[i] = 3*k+j with k is the k-th triangle of the ball and j the position of the point in the ball (0,1,2) */
int boulep(pMesh mesh, int start, int point , int** list)
{
	int i,j,compt=0;
	clock_t debut;
	clock_t fin;
	double difference;
	
	// clock just mesures the time to compare later the efficiency of this function with another
	debut= clock();
	// Parcour all triangles
	for(i=1 ; i <= mesh->nt ; i ++ )
	{
		// Parcour all 3 points of the triangle
		for ( j=0 ; j<=2 ; j++ )
		{
			if ( mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				compt++; // counts to define the malloc after
			}
		}
	}

	// define list	
	*list = (int*) malloc ( (compt+1) * sizeof(int)) ;
	compt = 0 ;
	for(i=1 ; i <= mesh->nt ; i ++ )
	{
		// Parcour all the 3 points of the triangle
		for ( j=0 ; j<=2 ; j++ )
		{
			if (  mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				// add number of triangle and number of vertice to list
				(*list)[compt] = 3*(i) + j ; 

				compt++;
			}
		}
	}
	fin = clock ();
	difference = difftime (fin, debut);
	//fprintf(stdout," temps d'éxécution %f \n",difference );
	return compt ; 
}


/* Function boule optimised!
   Searches the ball/triangles around the given point and saves them in list
   Does not walk anymore through all the triangles. Checks just in the 
   neighbourhood if the triangles belong to the point
   triangles are saved like this: list[i] = 3*k+j with k is the k-th triangle of the ball and j the position of the point in the ball (0,1,2) */
int boule_adj(pMesh mesh, int start, int point , int** list)
{
	int i,j,compt = 0;
	int triangle , triangle_av;
	*list = (int*)malloc(10 *sizeof(int)) ;
	int triangle_ori = start ;
	clock_t debut;
	clock_t fin;
	double difference;
	debut = clock ();
	do
	{
		for ( i = 0 ; i<3 ; i++ )
		{
			triangle =  mesh->adja[3*(triangle_ori-1)+1+i]/3 ; 
			if ( triangle != triangle_av )
			{
				for ( j=0 ; j<=2 ; j++ )
				{
					if (  mesh->tria[triangle].v[j] == mesh->tria[start].v[point] )
					{
						/* on remplit la list */
						(*list)[compt] = 3*triangle + j ;
						//fprintf(stdout," calcul %d \n",(*list)[compt]);
						triangle_av = triangle_ori ;
						triangle_ori = triangle ;
						compt++;
					}
				}
			}
		}
		
		
	
	}while ( triangle_ori != start  );
	fin = clock ();
	difference = difftime (fin, debut);
	//fprintf(stdout," temps d'éxécution %f \n",difference );
	return compt;
}



/* This function takes a array tab of the size of all points in the mesh
 For each point in the mesh it associates one possible triangle belonging to this point */
void hashTria(pMesh mesh, int *tab)
{
	int i, j, k;
	unsigned control = 0;
	//printf("ICI\n");
	for(i=1; i<=mesh->np; i++)
	{
		for(j=1; j<=mesh->nt && control == 0; j++)
		{
			for(k=0; k<3 && control == 0; k++)
			{
				if (mesh->point[i].c[0] == mesh->point[mesh->tria[j].v[k]].c[0] && mesh->point[i].c[1] == mesh->point[mesh->tria[j].v[k]].c[1] && mesh->point[i].c[2] == mesh->point[mesh->tria[j].v[k]].c[2])
				{
					tab[i] = j;
					control = 1;
					break;
				}
			}
		} 
		control = 0;
	}
}
