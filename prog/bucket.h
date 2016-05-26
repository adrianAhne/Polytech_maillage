#ifndef BUCKET_H
#define BUCKET_H


/* Structure bucket */
typedef struct {
  int size ;  // Nomber of subdivision for the axis  x,y,z 
  int *head ; // this is a tab of N³ in order to stock the points
  int *link ; // this tab has a size of mesh->np and this will contain the old_Points
} Bucket;
typedef Bucket* pBucket;

/* initialise the bucket*/
void init_bucket( pBucket bucket , pMesh mesh) ; 


/* Fill the 3 objects of a bucket
	First we will go through all the points of the mesh and calculate for each point:
						i = max(0,(int)(N*p->c[0])-1) 
						j = max(0,(int)(N*p->c[1])-1) 
						k = max(0,(int)(N*p->c[2])-1)
	Each of the integer is in the interval [0, N-1]
	Via this number we can associate each point with a key = (k*N+j)*N+i				
	With this key we go to bucket->head[key] and check : 
		- if  bucket->head[key] == 0, we put p in it 
		- if  bucket->head[key] > 0 , there is already a point (= Old_point) 
			So we put the old point at bucket->link[p] = bucket->head[key] 
			and we put the new point in head bucket->head[key] = p	   */
void fill_bucket( pBucket bucket , pMesh mesh ) ;

/* FUNCTION use_bucket_around  ( To search around ) 
		This function will search the neighbourhood for the point in the subdomains around the main subdomains.
		For this we will add one to the key or minus one
		Parameters : The bucket , the point , an increment
*/
int use_bucket_around(pBucket bucket,pPoint point,int increment, int* resultat,int key, int newkey);


/* 
	This function will use the coordinates of a point and associate the key.
	With this key, we will define the neighbourhood  of the point.
*/
int use_bucket( pBucket bucket , pMesh mesh ,  pPoint point , double increment )  ;


/* FUNCTION positive_boundingbox 
This function will translate the mesh using the 3D translation function in order to have a bounding box in the positive part of the axis (x,y,z)
*/
void positive_boundingbox( pMesh mesh , pPoint point);

/* 
	Determine the minimum point
*/
double point_min (pMesh mesh , char axis);

/* FUNCTION free_bucket */
void free_bucket (pBucket bucket) ;


/* Function uses bucket to return the value of the associated key to the case where the point is */
int bucket_retour_key( pBucket bucket , pMesh mesh ,  pPoint point , double increment );





#endif
