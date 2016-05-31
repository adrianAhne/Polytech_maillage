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




/* FUNCTION distbuck 
		This function will calculate the minimum distance between a point of the mesh A and the mesh B 
		PARAMETERS: a pointer to the point p of the mesh A, a pointer to the bucket of the mesh B, a pointer to the mesh B 
		RETURN : the distance d 
*/ 
double distbuck( pPoint p , pBucket bucket_meshB , pMesh meshB )
{
	/* we define d,dapp,iel,kel */
	double d = 10000.0, dapp = 10000.0 ;
	int kel = 0 , iel = 0 ;
	
	/* we save the box where the point p is in the bucket of the mesh B */
	int C = bucket_retour_key(  bucket_meshB , meshB , p , 0.0 );
	
	/* Now we search the nearest point p0 in the box of the bucket */
	/* We go throught all the points of the box and for each one we calculate the distance */
	Point p0 , pt = meshB->point[ bucket_meshB->head[C] ] ; ;
	double dist0 , d0 ;
	int cherche = bucket_meshB->head[C] ;
	/* First we calculate the distance for the first point */
	d0 = sqrt( ( pt.c[0] - p->c[0] ) * ( pt.c[0] - p->c[0] ) + ( pt.c[1] - p->c[1] ) * ( pt.c[1] - p->c[1] ) + ( pt.c[2] - p->c[2] ) * ( pt.c[2] - p->c[2] )) ;
	p0 = meshB->point[ bucket_meshB->head[C] ] ;
	p0.s =  cherche ;
	while(bucket_meshB->link[cherche])
	{
		/* we save the point */
		pt = meshB->point[ bucket_meshB->link[cherche] ] ;
		
		/* we calculate the distance between p0 and p */
		dist0 = ( pt.c[0] - p->c[0] ) * ( pt.c[0] - p->c[0] ) + ( pt.c[1] - p->c[1] ) * ( pt.c[1] - p->c[1] ) + ( pt.c[2] - p->c[2] ) * ( pt.c[2] - p->c[2] ) ;
		//printf ( " dist0 = %f  ",dist0) ;
		/* we save the smallest distance */
		if ( dist0 < d0 )
		{
			d0 = sqrt(dist0) ;
			p0 = pt ;
			p0.s =  bucket_meshB->link[cherche] ;
		}
		/* and we change the point */
		cherche = bucket_meshB->link[cherche] ;
		
	
	}
 /* printf ( " d0 = %f  ",d0) ;
	printf ( " p0.s = %d  ",p0.s) ;*/
	/* if d0 = 0 so d = 0 */
	if (!d0) 
		return d0 ; 
	/* Now we have the nearest point p0 of the bucket  */
	
	int test = 1 ;
	while ( test )
	{
		/* We make the list of the triangles around p0 */
		/* we search the first triangle containing the point p0 */
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
		/* We have the first triangle */
		/* Now we make the ball */
		int indice_point , indice_p0 ;
		double dk ;
		int** list_trianglep0 =(int**) malloc ( sizeof(int*)) ;
		int nb_trianglesp0 = boulep ( meshB, tr, indicept , list_trianglep0);

		/* For each triangle of the ball B0, we calculate the distance between the point p of the mesh A and the triangle of the ball B0 
			We save the smallest distance */
		for (i = 0 ; i < nb_trianglesp0 ; i++ )
		{
			/* we collect the index of the point p0 for each triangle k*/
			indice_point =  (*list_trianglep0)[i] % 3 ;
			/* we calculate the distance between the point p and each triangle k of the ball Bk */
			dk = distPointToTriangle(meshB, &meshB->tria[((*list_trianglep0)[i]-indice_point)/3], p);
		
			/* Now we save the minimum distance and the index of the triangle */
			if ( dk < d )
			{
				d = dk ;
				kel = ((*list_trianglep0)[i]-indice_point)/3 ;
				indice_p0 = indice_point;
			}
		}

		/* At this point, we have the triangle, around p0, close to the point p0 of the mesh A*/
		/* now for each point of this triangle different of p0 :
				we construct the ball around each point 
					we save the triangle closer to the point p 
						we compare the distance between this triangle and the point p and the distance between the triangle of the ball p0 and the point p 
						*/
		for( i = 0 ; i <=2 ; i++ )
		{		
			/* if different of p0 )*/		
			if ( i != indice_p0 )
			{
				/* we construct the ball around this new point: start is the triangle kel and the index is i */
				int indice_pointpk , j , indice_pkmin ;
				double d_pk ;
				int** list_trianglepk =(int**) malloc ( sizeof(int*)) ;
				int nb_trianglespk = boulep ( meshB, kel , i , list_trianglepk);
				
				/* At this point for each triangle of this ball we calculate the distance between the point of the mesh A and the triangles */
				for ( j = 0 ; j < nb_trianglespk ; j++ )
				{
					/* we collect the index for each triangle k*/
					indice_pointpk =  (*list_trianglepk)[j] % 3 ;				
					d_pk = distPointToTriangle(meshB, &meshB->tria[((*list_trianglepk)[j]-indice_pointpk)/3], p);
					
					/* we save the smallest distance and the index of the triangle closer */
					if( d_pk < dapp )
					{
						dapp = d_pk ;
						iel = ((*list_trianglepk)[j]-indice_pointpk)/3 ;
						indice_pkmin = indice_pointpk ;
					}
				}
				
				/* we compare the distances */
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
		
		/* free */
		free( (*list_trianglep0)) ;
		free(list_trianglep0) ;
	}
	
	return d ;
}





