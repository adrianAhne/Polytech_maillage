#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
//#include <omp.h> // parallel computing
#include "mesh.h"
#include "bucket.h"
#include "distance.h"
#include "ball.h"
#include "hash.h"


#define KA 7
#define KB 11
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


unsigned char idir[5]     = {0,1,2,0,1}; 


/*
Hedge struct doc:
na and nb are the endpoints of the edge

ia = min(na,nb)
ib = max(na,nb)

This objects Hedge are ordered in a global array

key = (KA*min(na,nb)+KB*max(na,nb))%hsize
KA = 7, KB = 11, hsize = mesh->np
*/



/* Construct the table hachage */
int hashHedge(pMesh mesh, Hedge *tab)
{
	pHedge ph;
	pTria  pt;
	int i,j,k,i1,i2,min,max;
  
	int hnxt = mesh->np + 1; // gives first available index in the array tab
	int hsize = mesh->np;
	int key;
	int adj,control=2;
  
  	// Calculate for each edge (na, nb) in the mesh the key key and the adjacent relation (3*i + j) adj by which this edge is seen in the triangle
	for(i=1; i<=mesh->nt; i++)
	{
    	pt = &mesh->tria[i];
    	if ( !pt->v[0] )  continue;
    	
    	for(j=0;j<3;j++)
    	{
    	    i1 = idir[j+1];
    		i2 = idir[j+2];
      		min = MIN(pt->v[i1], pt->v[i2]);
      		max = MAX(pt->v[i1], pt->v[i2]);
      
      		/* compute key */
	        key = (KA*min +KB*max)%hsize + 1;
      		
      		/* insert */
      		adj = 3*i+j;
      
      		// check is tab[key] exists
      		if (tab[key].adj1 == 0) 
      		{
      			tab[key].ia = min;
		        tab[key].ib = max;
        		tab[key].adj1 = adj;
      		}
      		// Search if one of the elements corresponds to (na, nb). If there is one, his data adj1 has already been filled.
      		// If not, go until the end of the chaine (until there is no more nxt), take the first available element in tab (which is the index of hnxt) and fill up the associated data Hedge.
			else { 
		        control =2 ;
		        do {
					if ( (tab[key].ia == min ) && (tab[key].ib == max)  )
        				{ tab[key].adj2 = adj; control =1 ;}
          
					key = tab[key].nxt;
        		}
		        while(( tab[key].nxt) && ( control==2 ) );
        
        
        		if(control==2) {
        			tab[hnxt].ia = min;
        			tab[hnxt].ib = max;
        			tab[hnxt].adj1 = adj;
        			tab[key].nxt = hnxt; // the last visited tab[key] points to the next element
        			hnxt++;
        		}
			}		
    	}
	}
  
	return 0;
  
}


/* Calculate neighbour/adjacent triangles of a given triangle (only at edges -> max 3 triangles) */
int setAdj(pMesh mesh, Hedge * tab)
{
  	pTria  pt;
  	Hedge  *ph;
	int i,j,i1,i2,min,max,adj;
	int key,control;
	int hsize = mesh->np; // gives first indice disponible in the tableau

	// 3 arrêt per triangle
	mesh->adja = (int *)calloc(3*mesh->nt+1,sizeof(int));

	// for chaque triangle, mesh->adja[3*(k-1)+1+j] will be the associated adjacent relation of the arrêt j of k, seen in the neighbour triangle
	for(i=1; i<=mesh->nt; i++)
	{
    	pt = &mesh->tria[i];
    	if ( !pt->v[0] )  continue;
    
		for(j=0; j<3; j++)
    	{
    		// Calculate the associated key to the edge (na, nb)
      		i1 = idir[j+1];
		    i2 = idir[j+2];
		    min = MIN(pt->v[i1], pt->v[i2]);
		    max = MAX(pt->v[i1], pt->v[i2]);
			key = (KA*min+KB*max)%hsize+1;
			
	   
	        // Parcour the chaine of the associated array tab at this key value, until what the object Hedge associated to (na, nb) is found
	        // One of the two data adj1, adj2 of this object is the adjacent relation of the edge in the triangle i. The other one is the searched relation to put in mesh->adja[3*(i-1) + 1 + j]
      		do{
      			if ( (tab[key].ia == min) &&  (tab[key].ib == max) )
      			{ 
      				adj = 3*i+j;
			        if( (tab[key].adj2) == adj)  
			        	mesh->adja[3*(i-1)+1+j] = tab[key].adj1;
			        else   
			        	mesh->adja[3*(i-1)+1+j] = tab[key].adj2;
          			break;
        		}
        		key = tab[key].nxt;
      		}
      		while(tab[key].nxt);
    	}
	}

	return 0;
}

/* Brute force method to find the index of a triangle for a given point */
int localiseTriangleBruteForce(pMesh mesh, pPoint point){
	int i,j;
	for(i=1; i<mesh->nt; i++){
		for(j=0; j<3; j++){
			if (point->c[0] == mesh->point[mesh->tria[i].v[j]].c[0])
				return i;
		}
	}
	return 0;
}

/* calculates coordinates barycentric u, v, 1-u-v for a given point p and saves them in the pointer cb */
void baryCoord(pMesh mesh, int triangle, pPoint p, double cb[3])
{
    Point a = mesh->point[ mesh->tria[triangle].v[0] ];
    Point b = mesh->point[ mesh->tria[triangle].v[1] ];
    Point c = mesh->point[ mesh->tria[triangle].v[2] ];

    Point v0, v1, v2;

    v0.c[0] = b.c[0] - a.c[0];
    v0.c[1] = b.c[1] - a.c[1];
    v0.c[2] = b.c[2] - a.c[2];

    v1.c[0] = c.c[0] - a.c[0];
    v1.c[1] = c.c[1] - a.c[1];
    v1.c[2] = c.c[2] - a.c[2];

    v2.c[0] = p->c[0] - a.c[0];
    v2.c[1] = p->c[1] - a.c[1];
    v2.c[2] = p->c[2] - a.c[2];

    double d00, d01, d11, d20, d21, denom;

    d00 = dotProduct3D(v0, v0);
    d01 = dotProduct3D(v0, v1);
    d11 = dotProduct3D(v1, v1);
    d20 = dotProduct3D(v2, v0);
    d21 = dotProduct3D(v2, v1);
    denom = d00 * d11 - d01 * d01;

    cb[0] = (d11 * d20 - d01 * d21) / denom; // u
    cb[1] = (d00 * d21 - d01 * d20) / denom; // v
    cb[2] = (1. - cb[0] - cb[1]); // 1 - u - v

}

/* find triangle who includes the point p by starting at a starttriangle and approaches the searched triangles via adjacents */
int locelt(pMesh mesh, int startTriangle, pPoint p, double cb[3])
{

	int triaRef = startTriangle;
	pTria t = &mesh->tria[triaRef];
	double distOld = distPointToTriangle(mesh, t, p);
	double distCurrent;
	int refClosest = 0;
	int i, it = 0;

	while(1)
	{
		// calculate barycentric coordinates and save them in array cb
		baryCoord(mesh, triaRef, p, cb);

		// if all coodinates positive, the current Triangle is returned
		if (cb[0] >= 0 && cb[1] >= 0 && cb[2] >= 0)
			return triaRef;
		else 
		{
			for(i=0; i<3; i++)
			{
				t = &mesh->tria[mesh->adja[3*(triaRef-1) + 1 + i]/3];
				assert(mesh->adja[3*(triaRef-1) + 1 + i]/3 <= 3*mesh->nt);
				//printf("Triangle numero %d\n", mesh->adja[3*(triaRef-1) + 1 + i]/3);
		
				// calculates average distance of point p to triangle t
				distCurrent = averageDistancePTT(mesh, t, p);
				
				if (distCurrent < distOld)
				{
					distOld = distCurrent;
					triaRef = mesh->adja[3*(triaRef-1) + 1 + i]/3;
					break;
				}
				if (i == 2)
				{
					// this means that unter the neighbours of a triangle, we did not find a neighbour who is the closest therefore we return with -2
					return -2; 
				}
			}
		}
		
		it++;

		// prevents only that the program does not stock in an endless loop
		if (it > 100000)
			return -1;
	}
}


/*	 Following the algorithm given in github/doc/main.pdf
	 Compute the distance between a point and a surface using buckets and adjacences
	 VertToTria gives an array which tells us to which triangle a given point belongs to */
double distanceUsingBucket(pMesh mesh, pPoint p, int *VertToTria , Bucket* bucket )
{
 	unsigned N = 0;
	double d = 100000000.0, dapp = 10000.0, d0, dk, d_pk;
	int kel = 0, iel = 0, C, p0, pk, i,j,k, start, pointN, ind ;
	int TriaInBoule, TriaInBoulePK;


	int indice_point ;
	int real_indic ;

	//positive_boundingbox( mesh , &point );

	
	// Find the cell grid C to which point belongs 
	// returns number of bucket
	C = bucket_retour_key( bucket , mesh, p , 0.0 ); // doesn't work
	printf("Cell grid of the point : %d\n", C);

	// TO DO: how to choose p0???
	// check for each point in the bucket the distance to the given point and take the minimum as p0 (Norm2)
	p0 = bucket->head[C]; // false : one point in the cell grid with it's distance to point
	//ALEXANDRE: Comment peux-je avoir accès aux elements du bucket? 
	// Réponse: je vais te créer un tableau avec tout les points d'une case 
	int* points ; // le tableau d'indice des points

	int cherche = bucket->head[C] ;
	int compte = 0 ;
	
	while ( bucket->link[cherche] ) 
	{
		cherche = bucket->link[cherche] ;
		compte ++ ;
	}
	points =(int*) malloc((compte+1) * sizeof(int)	 ) ;
	/* fill up the array of indice of points */
	points[0] = bucket->head[C] ;
	cherche = bucket->head[C] ;
	for( i = 1 ; i <= compte ; i++ ) 
	{
		points[i] = bucket->link[cherche] ;
		//fprintf(stdout,"  Points[%d]  x = %f  y = %f  z =  %f \n",i, mesh->point[points[i]].c[0],mesh->point[points[i]].c[1],mesh->point[points[i]].c[2]);
		cherche  = bucket->link[cherche] ;
	} 
	
	printf("compte = %d \n" , compte);
	/* ICI tu as le tableau remplis tu peux calculer la distance : bon courage */

	
	// Starting from C explore the bucket and find the vertex triangulation p0 closer to p and retain the distance d0 = d(p,p0)
	double dist0 = 100000000000000;
	double dCurrent;
	for(i=0; i <= compte; i++)
	{
		// 3-dimensional case: calculate distance between point in bucket and given point p via the norm2 : dCurrent = sqrt(x_bucket - x_point)^2 + (y_bucket - y_point)^2 + (z_bucket - z_point)^2
		dCurrent = sqrt((mesh->point[points[i]].c[0] - p->c[0])*(mesh->point[points[i]].c[0] - p->c[0]) + (mesh->point[points[i]].c[1] - p->c[1])*(mesh->point[points[i]].c[1] - p->c[1]) + (mesh->point[points[i]].c[2] - p->c[2])*(mesh->point[points[i]].c[2] - p->c[2]));

		if (dCurrent < dist0)
		{
			dist0 = dCurrent;
			p0 = points[i]; // integer of current point in bucket	
		}
	}

	d0 = dist0;
		printf("d0 = %f \n",d0 ) ;

	if ( !d0 )
		return d0 ;
	
	// startin from C explore the bucket and find the vertex triangulation p0 closer to p and retain the distance d0 = d(p,p0)
	int thepoint = p0 ;
	printf(" p0 = %d\n",p0);
	while(p0)
	{
		int** list = (int**)malloc(sizeof(int*)) ; // list of triangles in the ball B(p0)
	

		// get list of all triangles in the ball B(p0) of p0
		TriaInBoule = boulep(mesh, VertToTria[p0], p0 , list);

		// for each triangle K_k' in the ball B(p0)
		for(k=0; k < TriaInBoule; k++)
		{
				
			//printf(" %d \n",k);
			//printf(" *(list)[k]/3 = %d \n",(*list)[k]/3 ) ; 
			// calculate distance d_k = d(p, K_k');
			int indice_point =  (*list)[k] % 3;
			//printf("Liste des triangles autour :  \n " );
			/*for (i=0;i<TriaInBoule;i++)
				printf("triangle %d = %d \n ", i,((*list)[i] - indice_point)/3 );
			*/

			dk = distPointToTriangle(mesh, &mesh->tria[((*list)[k]-indice_point)/3], p);
			//printf(" dk = %f \n",dk);

			if (dk < d0)
			{
				//update the distance d= dk
				d0 = dk;
				// store index of current triangle kel = k
				kel = ((*list)[k]-indice_point)/3 ;
				real_indic = indice_point ;			
			}

		}
		free(list);
			
		
		// for each vertex pk belonging to the triangle kel
		for(i=0; i <= 2; i++)
		{

			//printf("p0 = %d et mesh->tria[kel].v[i] = %d \n" , p0 , mesh->tria[kel].v[i] );
			// check that current point is not the point, which was already handled in the for loop before. Consider only the two vertices of the triangle which are left
			if (mesh->tria[kel].v[i] != p0)
			{
				int** listLocal = (int**)malloc(sizeof(int*)); // list of triangles in the ball B(pk)
				//printf("ICI\n") ; 

				pk = mesh->tria[kel].v[i] ;

				// get list of all triangles in the ball B(pk) of pk
				TriaInBoulePK = boulep(mesh, kel, i, listLocal);
				/*for (i=0;i<TriaInBoulePK && kel!= (* listLocal)[i] ;i++)
					printf("triangle %d = %d \n ", i,(* listLocal)[i] );*/


				// for each triangle K_pk' in the ball B(pk)
		
				for (j=0; j < TriaInBoulePK  ;j++)
				{
						
					int indice_point = (*listLocal)[j] %3;

					if ( (*listLocal)[j] != kel ) 
					{
						//compute the distance d_pk = d(p, K_pk')
						d_pk = distPointToTriangle(mesh, &mesh->tria[((*listLocal)[j]-indice_point)/3], p);

			

						if(d_pk < dapp)
						{
							// update the distance
							dapp = d_pk;
							
							// store the index of current triangle iel = j
							iel = j;
						}
					}
				}
					free(listLocal);
				//printf("dapp = %f  d = %f\n ", dapp,d) ;

		
			if (dapp < d0)
			{
				p0 = pk; 
				if (p0 == thepoint )
					p0 = 0;
			}
			else
				p0 = 0;
		
			
			}
			
		}	

	}
	



	return d0;
}


double distbuck( pPoint p , pBucket bucket_meshB , pMesh meshB )
{
	/* on définit d,dapp,iel,kel */
	double d = 10000.0, dapp = 10000.0 ;
	int kel = 0 , iel = 0 ;
	
	/* On stocke la case ou le point p est dans le bucket du maillage B */
	int C = bucket_retour_key(  bucket_meshB , meshB , p , 0.0 );
	
	/* Maintenant on cherche le point p0 de la case du bucket le plus proche de p */
	/* On parcours tous les points de la case et pour chacun d'entre eux on calcule la distance */
	Point p0 , pt = meshB->point[ bucket_meshB->head[C] ] ; ;
	double dist0 , d0 ;
	int cherche = bucket_meshB->head[C] ;
	/* on calcule la distance pour le premier point */
	d0 = sqrt( ( pt.c[0] - p->c[0] ) * ( pt.c[0] - p->c[0] ) + ( pt.c[1] - p->c[1] ) * ( pt.c[1] - p->c[1] ) + ( pt.c[2] - p->c[2] ) * ( pt.c[2] - p->c[2] )) ;
	p0 = meshB->point[ bucket_meshB->head[C] ] ;
	p0.s =  cherche ;
	while(bucket_meshB->link[cherche])
	{
		/* on récupère le point */
		pt = meshB->point[ bucket_meshB->link[cherche] ] ;
		
		/* on calcule la distance entre p0 et p */
		dist0 = ( pt.c[0] - p->c[0] ) * ( pt.c[0] - p->c[0] ) + ( pt.c[1] - p->c[1] ) * ( pt.c[1] - p->c[1] ) + ( pt.c[2] - p->c[2] ) * ( pt.c[2] - p->c[2] ) ;
		//printf ( " dist0 = %f  ",dist0) ;
		/* on stocke la plus petite */
		if ( dist0 < d0 )
		{
			d0 = sqrt(dist0) ;
			p0 = pt ;
			p0.s =  bucket_meshB->link[cherche] ;
		}
		/* on change de point */
		cherche = bucket_meshB->link[cherche] ;
		
	
	}
 /* printf ( " d0 = %f  ",d0) ;
	printf ( " p0.s = %d  ",p0.s) ;*/
	if (!d0) 
		return d0 ; 
	/* A présent on a le point p0 de la case le plus proche de p */
	/* On rentre dans la boucle */
	int test = 1 ;
	while ( test )
	{
		/* On fait la liste des triangles autour de p0 */
		/* Tout d'abord on cherche le premier triangle qui contient p0 pour nous servir de départ à la boule */
		int tr = 0 ,i = 1 , indicept;
		while ( !tr )
		{
			if( meshB->tria[i].v[0] == p0.s )
			{
				tr = i ;
			 	indicept = 0 ;
			}
			if( meshB->tria[i].v[1] == p0.s )
			{
				tr = i ;
				indicept = 1 ;
			}
			if( meshB->tria[i].v[2] == p0.s )
			{
				tr = i ;
				indicept = 2 ;
			}
			i++;
		}
		/* on a donc le premier triangle et l'indice du point correspondant */
		/* A présent on peut faire la boule */
		int indice_point , indice_p0 ;
		double dk ;
		int** list_trianglep0 =(int**) malloc ( sizeof(int*)) ;
		int nb_trianglesp0 = boulep ( meshB, tr, indicept , list_trianglep0);
		/*printf("Liste des triangles autour :  \n " );
		for (i=0;i < nb_trianglesp0;i++)
		{
			indice_point =  (*list_trianglep0)[i] % 3;
			printf("triangle %d = %d \n ", i,((*list_trianglep0)[i] - indice_point)/3 );
		}	*/
		/* A présent pour chaque triangle de la boule B0 on calcule la distance entre le point p du maillage A et les triangles de la boule autour de p0 
			On stocke la plus petite distance */
		for (i = 0 ; i < nb_trianglesp0 ; i++ )
		{
			/* on récupère l'indice du point p0 pour chaque triangle k*/
			indice_point =  (*list_trianglep0)[i] % 3 ;
			//printf ( " indice_point = %d ", indice_point ) ;
			/* on calcule à présent la distance entre le point p et chaque triangle k*/
			dk = distPointToTriangle(meshB, &meshB->tria[((*list_trianglep0)[i]-indice_point)/3], p);
		
			/* on stocke à présent la distance minimale et l'indice du triangle correspondant */
			if ( dk < d )
			{
				d = dk ;
				kel = ((*list_trianglep0)[i]-indice_point)/3 ;
				indice_p0 = indice_point;
			}
		}
	//	printf( " d = %f  ",d); 
		/* A présent on a le triangle, autour du point p0, le plus proche du point du meshA  */
		/* Maintenant pour tous les points de ce triangles différents de p0 :
				on construit la boule autour de ces points 
					on stocke le triangle le plus proche de p 
						on compare la distance entre le point et ce triangle et celle entre le point et le triangle de la boule autour de p0 */
		for( i = 0 ; i <=2 ; i++ )
		{		
			/* Si différent de p0 )*/		
			if ( i != indice_p0 )
			{
				/* on construit la boule autour de ce nouveau point : le triangle de départ et le triangle kel et l'indice du point est i*/
				int indice_pointpk , j , indice_pkmin ;
				double d_pk ;
				int** list_trianglepk =(int**) malloc ( sizeof(int*)) ;
				int nb_trianglespk = boulep ( meshB, kel , i , list_trianglepk);
				
				/* A present pour chaque triangle de cette boule on calcule la distance entre le point du maillage A et les triangles */
				for ( j = 0 ; j < nb_trianglespk ; j++ )
				{
					/* on récupère l'indice du point pour chaque triangle k*/
					indice_pointpk =  (*list_trianglepk)[j] % 3 ;				
					d_pk = distPointToTriangle(meshB, &meshB->tria[((*list_trianglepk)[j]-indice_pointpk)/3], p);
					
					/* on stocke la plus petite distance et l'indice du triangle le plus proche */
					if( d_pk < dapp )
					{
						dapp = d_pk ;
						iel = ((*list_trianglepk)[j]-indice_pointpk)/3 ;
						indice_pkmin = indice_pointpk ;
					}
				}
				
				/* On compare à présent la distance obtenu par le triangle de la boule Bk à celui de la boule B0 */
				if ( dapp < d ) 
				{
					d = dapp ;
					test = meshB->tria[iel].v[indice_pkmin] ; 
					p0 = meshB->point[test] ;
					p0.s = test ;
				}
				else
					test = 0 ;
				free( (*list_trianglepk)) ;
				free(list_trianglepk) ;
				
				
			}	
	
		}
		
		/* on libère les listes avant de recommencer */
		free( (*list_trianglep0)) ;
		free(list_trianglep0) ;
	
	}
//	printf("d return = %f \n",d);
	return d ;
	
}





