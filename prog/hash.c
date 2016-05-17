#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // you need to add -lm for link edition
//#include <omp.h> // parallel computing
#include "mesh.h"
#include "hash.h"
#include "distanceMeshFunctions.h"
#include "bucket.h"
#include "boule.h"

#define KA 7
#define KB 11
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


unsigned char idir[5]     = {0,1,2,0,1}; 

int hashHedge(pMesh mesh, Hedge *tab)
{
  pHedge ph;
  pTria  pt;
  int i,j,k,i1,i2,min,max;
  
  int hnxt = mesh->np + 1;
  int hsize = mesh->np;
  int key;
  int adj,control=2;
  
  
  for(i=1; i<=mesh->nt;i++)
  {
    pt = &mesh->tria[i];
    if ( !pt->v[0] )  continue;
    
    for(j=0;j<3;j++)
    {       i1 = idir[j+1];
      i2 = idir[j+2];
      min = MIN(pt->v[i1], pt->v[i2]);
      max = MAX(pt->v[i1], pt->v[i2]);
      /* compute key */
      key = (KA*min +KB*max)%hsize + 1;
      /* insert */
      adj = 3*i+j;
      
      
      if (tab[key].adj1 == 0) // si tab[key] est visité pour la premiere fois
      {
        tab[key].ia = min;
        tab[key].ib = max;
        tab[key].adj1 = adj;
      }
      
      else { // tab[key].adj1 a déjà été rempli
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
          tab[key].nxt = hnxt;
          hnxt++;
        }
      }
    }
  }
  
  
//  /* test debugging */
//  /* On cherche les trois voisins du triangle 9629 */
//  pt = &mesh->tria[9629];
//  for(j=0;j<3;j++)
//  {  i1 = idir[j+1];
//    i2 = idir[j+2];
//    fprintf(stdout,"  -- j %d  i1 %d  i2 %d \n",j,i1,i2);
//    min = MIN(pt->v[i1], pt->v[i2]);
//    max = MAX(pt->v[i1], pt->v[i2]);
//    /* compute key */
//    key = (KA*min +KB*max)%hsize + 1;
//    /* insert */
//   
//     
//    do {
//      if ( (tab[key].ia == min) &&  (tab[key].ib == max) )
//
//    fprintf(stdout,"  -- ad1 %d  adj2 %d  \n",tab[key].adj1/3,tab[key].adj2/3);
//
//       key = tab[key].nxt;
//        
//    }
//    
//    while(tab[key].nxt);
//  
//  }
  
  
  return 0;
  
}


// Calculate neighbour triangles of a given triangle (only at edges -> max 3 triangles)
int setAdj(pMesh mesh, Hedge * tab)
{
  pTria  pt;
  Hedge  *ph;
	int i,j,i1,i2,min,max,adj;
	int key,control;
	int hsize = mesh->np; // gives first indice disponible in the tableau

	// 3 arrêt per triangle
	mesh->adja = (int *)calloc(3*mesh->nt+1,sizeof(int));

	for(i=1; i<=mesh->nt; i++)
	{
    pt = &mesh->tria[i];
    if ( !pt->v[0] )  continue;
    
		for(j=0; j<3; j++)

    {
      i1 = idir[j+1];
      i2 = idir[j+2];
      min = MIN(pt->v[i1], pt->v[i2]);
      max = MAX(pt->v[i1], pt->v[i2]);
			key = (KA*min+KB*max)%hsize+1;
			
      // On parcourt la chaine de tab jusqu'à ce que l'objet correspondant à l'arrete soit trouvé
      do{
       if ( (tab[key].ia == min) &&  (tab[key].ib == max) )
        { adj = 3*i+j;
          if( (tab[key].adj2) == adj)  mesh->adja[3*(i-1)+1+j] = tab[key].adj1;
          else   mesh->adja[3*(i-1)+1+j] = tab[key].adj2;
          break;
        }
        key = tab[key].nxt;
      }
      while(tab[key].nxt);
    }
  }
  
  /* Test debugging: on cherche les trois voisins du triangle 34 
  int iadr,iel, *adja;
  iadr = (34-1)*3 + 1;
  adja = &mesh->adja[iadr];
  
  iel = (adja[0]) / 3 ;
  printf("Tria = %d\n", iel);
  iel = (adja[1]) / 3 ;
  printf("Tria = %d\n", iel);
  iel = (adja[2] ) / 3  ;
  printf("Tria = %d\n", iel);
  */
  return 0;

}


int localiseTriangleBruteForce(pMesh mesh, pPoint point){
	int i,j;
	for(i=1; i<mesh->nt; i++){
		for(j=0; j<3; j++){
			if (point->c[0] == mesh->point[mesh->tria[i].v[j]].c[0])
			{
				return i;
			}
		}
	}
	return 0;
}

// calculates coordinates barycentric u, v, 1-u-v
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

// find triangle who includes the point p
int locelt(pMesh mesh, int startTriangle, pPoint p, double cb[3])
{
	// distPointToTriangle(pMesh mesh, pTria tria, pPoint P0)
	int triaRef = startTriangle;
	pTria t = &mesh->tria[triaRef];
	printf("Start point = %d\n", t->v[0]);
	double distOld = distPointToTriangle(mesh, t, p);
	double distCurrent;
	int refClosest = 0;
	int i, it = 0;

	while(1)
	{
		baryCoord(mesh, triaRef, p, cb);
		printf("triaref : %d, cb : %f %f %f\n", triaRef, cb[0], cb[1], cb[2]);
		if (cb[0] >= 0 && cb[1] >= 0 && cb[2] >= 0)
		{
			printf("%d\n", triaRef);
			return triaRef;
		} else {
			for(i=0; i<3; i++)
			{

				t = &mesh->tria[mesh->adja[3*(triaRef-1) + 1 + i]/3];
				assert(mesh->adja[3*(triaRef-1) + 1 + i]/3 <= 3*mesh->nt);
				printf("Triangle numero %d\n", mesh->adja[3*(triaRef-1) + 1 + i]/3);
				//distCurrent = distPointToTriangle(mesh, t, p);
				distCurrent = averageDistancePTT(mesh, t, p);
				
				if (distCurrent < distOld)
				{
					distOld = distCurrent;
					triaRef = mesh->adja[3*(triaRef-1) + 1 + i]/3;
					printf("dist Old : %f triaRef : %d\n", distOld, triaRef);
					break;
				}
				if (i == 2)
				{
					return -2; 
					// c'est à dire que parmis les voisins d'un triangle, il n'a pas trouvé un voisin qui est plus proche
					// donc s'il continue il va s'éloigner du traingle à chercher
				}
			}
		}
		it++;
		if (it > 100000)
		{
			return -1;
		}
	}
}


// Following the algorithm given in github/doc/main.pdf
// Compute the distance between a point and a surface using buckets and adjacences
double distanceUsingBucket(pMesh mesh, pPoint p)
{
 	unsigned N = 0;
	double d = 10000.0, dapp = 10000.0, d0, dk, d_pk;
	int kel = 0, iel = 0, C, p0, pk, i,j,k, start, pointN ;
	int TriaInBoule, TriaInBoulePK;
	int** list = (int**)malloc(sizeof(int*)) ; // list of triangles in the ball B(p0)
	int** listLocal = (int**)malloc(sizeof(int*)); // list of triangles in the ball B(pk)

	// Hash et adja relations
	Hedge *tab = (Hedge*)calloc(3*mesh->nt+1,sizeof(Hedge));
	hashHedge(mesh, tab);
	setAdj(mesh, tab);
	
	int *hashT = (int*)calloc(mesh->np+1, sizeof(int));
	hashTria(mesh, hashT);

	// Bucket creation
	Bucket bucket;
	Point point;
	// point aleatoire
	point.c[0] =  .3;
	point.c[1] =  .4;
	point.c[2] =  .5;
	positive_boundingbox( mesh , &point );
	fprintf(stdout,"\n  -- Creation bucket MESH \n\nPlease type the number of subdivision : \n");
	fflush(stdin);	// just to be able to do a scanf without conflicts
	fscanf(stdin,"%d",&N);
	bucket.size = N ;
	init_bucket( &bucket , mesh); 
	fill_bucket( &bucket , mesh );

	// Find the cell grid C to which point belongs 
	// returns number of bucket
	C = bucket_retour_key( &bucket , mesh, &point , 0.0 ); // doesn't work
	printf("Cell grid of the point : %d\n", C);

	// TO DO: how to choose p0???
	p0 = bucket.head[C]; // false : one point in the cell grid with it's distance to point
	d0 = 1000.0; // false
	d0 = 
	
	// startin from C explore the bucket and find the vertex triangulation p0 closer to p and retain the distance d0 = d(p,p0)
	while(p0)
	{

		// get list of all triangles in the ball B(p0) of p0
		// TO DO: give startTriangle and startpoint
		TriaInBoule = boulep(mesh, start, pointN , list);
		
		// for each triangle K_k' in the ball B(p0)
		for(k=0; k < TriaInBoule; k++)
		{
			// calculate distance d_k = d(p, K_k');
			dk = distPointToTriangle(mesh, *list[k]/3, p0);
			if (dk < d)
			{
				//update the distance d= dk
				d = dk;
				// store index of current triangle kel = k
				kel = k;
			}
		}
		
		// for each vertex pk belonging to the trianlge kel
		for(i=0; i < 2; i++)
		{
			// get list of all triangles in the ball B(pk) of pk
			// TO DO: give startTriangle
			TriaInBoulePK = boulep(mesh, start, mesh->tria[kel].v[i], listLocal);
			// for each triangle K_pk' in the ball B(pk) 
			for (j=0; j < TriaInBoulePK; j++)
			{
				//compute the distance d_pk = d(p, K_pk')
				d_pk = distPointToTriangle(mesh, *list[j]/3, p);
				
				if(d_pk < dapp)
				{
					// update the distance
					dapp = d_pk;
					
					// store the index of current triangle iel = j
					iel = j;
				}
			}
			
			if (dapp > d)
				p0 = pk;
			else
				p0 = 0;
		}
	}
	

	free_bucket (&bucket);
	free(list);
	free(listLocal);

	return d;
}





