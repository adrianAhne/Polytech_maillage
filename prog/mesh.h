

#define NPMAX 10000

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


typedef struct {
  double    c[3],ref;
} Point;
typedef Point * pPoint;

typedef struct {
  int       v[2],ref,tag;
} Edge;
typedef Edge * pEdge;
  
typedef struct {
  int       v[3],ref;
} Tria;
typedef Tria * pTria;


// Adrian et Tupac: structure for edge normale
typedef struct {
  pPoint    point[3];
  pEdge		norm[3];
  int		ref;
} TriaNorm;
typedef TriaNorm * pTriaNorm;


typedef struct {
  int       np,na,nt,nr,ver,dim,mark;
  char     *namein,*nameout;
	double   *sol,o[3],rad;
  pPoint    point;
  pEdge     edge;
  pTria     tria;
  pTriaNorm triaNorm;
} Mesh;
typedef Mesh * pMesh;






/* read mesh */
int loadMesh(pMesh mesh);
 
/** Save mesh data */
int saveMesh(pMesh mesh);

/* prototypes */
int saveSol(pMesh mesh,int it);



 

