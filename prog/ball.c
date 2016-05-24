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




int boulep(pMesh mesh, int start, int point , int** list)
{
	int i,j,compt=0;
	clock_t debut;
	clock_t fin;
	double difference;
	
	//fprintf(stdout," Le point considéré est le %d \n",mesh->tria[start].v[point]);
	debut= clock();
	/* Méthode actuelle méthode naïve avec parcours de tout les triangles*/
	for(i=1 ; i <= mesh->nt ; i ++ )
	{
		/* parcours des 3 points de chaques triangles */
		for ( j=0 ; j<=2 ; j++ )
		{
			if ( mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				//fprintf(stdout," face %d \n",i);
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
		for ( j=0 ; j<=2 ; j++ )
		{
			if (  mesh->tria[i].v[j] == mesh->tria[start].v[point] )
			{
				/* on remplit la list */
			//	printf( "i  %d \n",i);
				(*list)[compt] = 3*(i) + j ;

				//fprintf(stdout,"*(list[%d] = 3*%d + %d = %d \n",compt, i, j,(*list)[compt]);
				compt++;
			}
		}
	}
	fin = clock ();
	difference = difftime (fin, debut);
	//fprintf(stdout," temps d'éxécution %f \n",difference );
	return compt ; 
	
	
}



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
			//fprintf(stdout," ici  " );
			triangle =  mesh->adja[3*(triangle_ori-1)+1+i]/3 ; 
			//fprintf(stdout," triangle_ori = %d " , triangle_ori ) ;
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



// This function takes a array of all the points in the mesh
// It completes the array it in order to associate one triangle to each point
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
