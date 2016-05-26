#ifndef BUCKET_H
#define BUCKET_H


/* Structure bucket */
typedef struct {
  int size ;  // Nomber of subdivision for the axis  x,y,z 
  int *head ; // this is a tab of N³ in order to stock the points
  int *link ; // this tab has a size of mesh->np and this will contain the old_Points
} Bucket;
typedef Bucket* pBucket;

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
void fill_bucket( pBucket bucket , pMesh mesh ) ;

/* FUNCTION use_bucket_around  ( To search around ) 
		This function will search the neighbourhood for the point in the subdomains around the main subdomains.
		For this we will add one to the key or minus one
		Parameters : The bucket , the point , an increment
*/
int use_bucket_around(pBucket bucket,pPoint point,int increment, int* resultat,int key, int newkey);

/* FUNCTION use_bucket 
		This function will use the coordinates of a point and associate the key.
		With this key, we will define the neighbourhood  of the point.
		Parameters : the bucket and a point 
		Return : a tab of the nearest points 
*/

int use_bucket( pBucket bucket , pMesh mesh ,  pPoint point , double increment )  ;


/* FUNCTION positive_boundingbox 
This function will translate the mesh using the 3D translation function in order to have a bounding box in the positive part of the axis (x,y,z)
PARAMETERS : a pointer to a mesh 
*/
void positive_boundingbox( pMesh mesh , pPoint point);

/* FUNCTION point_min 
In order to determine the minimum point  
PARAMETERS : a pMesh , a char for the axis
*/
double point_min (pMesh mesh , char axis);

/* FUNCTION free_bucket */
void free_bucket (pBucket bucket) ;

/* Fonction ou la fonction use bucket renvoi la valeur de la clé associé à la case ou se situe le point */
int bucket_retour_key( pBucket bucket , pMesh mesh ,  pPoint point , double increment );





#endif
