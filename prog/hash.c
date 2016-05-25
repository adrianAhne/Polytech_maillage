#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
//#include <omp.h> // parallel computing
#include "mesh.h"
#include "bucket.h"
#include "hash.h"
#include "distance.h"
#include "ball.h"

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
	double d = 10000.0, dapp = 10000.0, d0, dk, d_pk;
	int kel = 0, iel = 0, C, p0, pk, i,j,k, start, pointN, ind ;
	int TriaInBoule, TriaInBoulePK;
	int** list = (int**)malloc(sizeof(int*)) ; // list of triangles in the ball B(p0)
	int** listLocal = (int**)malloc(sizeof(int*)); // list of triangles in the ball B(pk)
	int indice_point ;
	int real_indic ;

	//positive_boundingbox( mesh , &point );

	
	// Find the cell grid C to which point belongs 
	// returns number of bucket
	C = bucket_retour_key( bucket , mesh, p , 0.0 );
	

	//p0 = bucket->head[C]; // false : one point in the cell grid with it's distance to point


	// saves all the points of the bucket in the array points
	int* points ; 
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
		cherche  = bucket->link[cherche] ;
	} 
	
	
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
	
	
	while(p0)
	{

		// get list of all triangles in the ball B(p0) of p0
		TriaInBoule = boulep(mesh, VertToTria[p0], p0 , list);

		// for each triangle K_k' in the ball B(p0)
		for(k=0; k < TriaInBoule; k++)
		{

			// calculate distance d_k = d(p, K_k');
			int indice_point =  (*list)[k] % 3;

			dk = distPointToTriangle(mesh, &mesh->tria[((*list)[k]-indice_point)/3], p);
			

			if (dk < d)
			{
				//update the distance d= dk
				d = dk;
				// store index of current triangle kel = k
				kel = ((*list)[k]-indice_point)/3 ;
				real_indic = indice_point ;			
			}
		}
			
		
		// for each vertex pk belonging to the triangle kel
		for(i=0; i <= 2; i++)
		{
			//printf("p0 = %d et mesh->tria[kel].v[i] = %d \n" , p0 , mesh->tria[kel].v[i] );
			// check that current point is not the point, which was already handled in the for loop before. Consider only the two vertices of the triangle which are left
			if (p0 != mesh->tria[kel].v[i])
			{
	
				pk = mesh->tria[kel].v[i] ;

				// get list of all triangles in the ball B(pk) of pk
				TriaInBoulePK = boule_adj(mesh, kel, i, listLocal);

				// for each triangle K_pk' in the ball B(pk)
				for (j=0; j < TriaInBoulePK && kel != (*listLocal)[j] ; j++)
				{
						int indice_point = (*listLocal)[j] %3;

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
			
				//printf("dapp = %f  d = %f\n ", dapp,d) ;
				
				if (dapp < d)
					p0 = pk; 
				else
					p0 = 0;
			}	
		}
			
			
		
	}
	
	free(list);
	free(listLocal);

	return d;
}






