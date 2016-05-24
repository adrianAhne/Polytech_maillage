#ifndef BALL_H
#define BALL_H




/* FONCTION boulep 

Cette fonction permet de retourner la boule autour d'un point. 
Cette boule est constituée de l'ensemble des points partageant ce point/sommet

Parametres : 
	- mesh  : le pointeur vers le maillage
	- start : l'indice du triangle 
	- point : le point du triangle start que l'on veut comme centre de la boule
	- list  : la liste des indices des triangles autour du sommet point
	
	return : le nombre de triangle ayant le point comme sommet ;
La liste des triangles est classée de la manière suivante: liste[i] = 3*k+j avec k le ième triangle de la boule et j la position du point dans la boule (0,1,2)
*/

int boulep(pMesh mesh, int start, int point , int** list);


int boule_adj(pMesh mesh, int start, int point , int** list);

void hashTria(pMesh mesh, int *tab);

#endif

