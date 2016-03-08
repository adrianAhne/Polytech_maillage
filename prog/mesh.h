

#define NPMAX 10000

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


typedef struct {
  double    c[3],n[3],ref;
} Point;
typedef Point * pPoint;

typedef struct {
  double    n[3];
} Normal;
typedef Normal * pNormal; 

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
  // Normal vertice defined by 2 points : norm start and norm end
  double  n[3];
  double  weight;
} TriaNorm;
typedef TriaNorm * pTriaNorm;


typedef struct {
  int       np,na,nt,nr,nn,ver,dim,mark;
  char     *namein,*nameout;
	double   *sol,o[3],rad;
  pPoint    point;
  pNormal   Normal;
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



 

