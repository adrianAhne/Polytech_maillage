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

#include "boule.h"



int boulep(pMesh mesh, int start, int point , int** list)
{
	int i,j,compt=0;
	fprintf(stdout," Le point considéré est le %d \n",mesh->tria[start].v[point]);
	/* Méthode actuelle méthode naïve avec parcours de tout les triangles*/
	for(i=1 ; i <= mesh->nt ; i ++ )
	{
		/* parcours des 3 points de chaques triangles */
		for ( j=0 ; j<2 ; j++ )
		{
			if ( mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				fprintf(stdout," face %d \n",i);
				/* on incremente le compteur qui nous servira pour le malloc du tableau list*/
				compt++;
			}
		}
	}
	
	*list = (int*) malloc ( (compt+1) * sizeof(int)) ;
	compt = 0 ;
	for(i=1 ; i <= mesh->nt ; i ++ )
	{
		/* parcours des 3 points de chaques triangles */
		for ( j=0 ; j<2 ; j++ )
		{
			if (  mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				/* on remplit la list */
				(*list)[compt] = 3*(i) + j ;
				fprintf(stdout," calcul %d \n",(*list)[compt]);
				compt++;
			}
		}
	}
	return compt ; 
}
