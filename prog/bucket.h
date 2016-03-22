#ifndef BUCKET_H
#define BUCKET_H


/* Structure bucket */
typedef struct {
  int size ;  // Nomber of subdivision for the axis  x,y,z 
  int *head ; // this is a tab of N³ in order to stock the points
  int *link ; // this tab has a size of mesh->np and this will contain the old_Points
} Bucket;
typedef Bucket * pBucket;

/* FUNCTION init_bucket 
			Parameters : a pointer for a bucket 
			Goal : Initialize head and link 
*/
void init_bucket( pBucket bucket , pMesh mesh) ; 

/* FUNCTION fill_bucket
			Parameters : a pointer for the bucket , one for the mesh , the number of subdivision  
			Goal : fill the 3 objects of a bucket 
			
			First we will go trough all the points of the mesh and calculate i,j,k for each point p
				
						i = max(0,(int)(N*p->c[0])-1) 
						j = max(0,(int)(N*p->c[1])-1) 
						k = max(0,(int)(N*p->c[2])-1) 
				Each of the int is in the interval [0;N-1]
				With this numbers we can associate each point with a key
			
				key = (k*N+j)*N+i
			
				With this key we go at bucket->head[key] : 
						- if  bucket->head[key] == 0, we put p in it 
						- if  bucket->head[key] > 0 , there is already a point (= Old_point) 
									So we put the old point at bucket->link[p] = bucket->head[key] 
									and we put the new point in head bucket->head[key] = p
											
			End boucle
			
*/
			
void fill_bucket( pBucket bucket , pMesh mesh , int N) ;

/* FUNCTION free_bucket */
void free_bucket (pBucket bucket) ;







#endif
