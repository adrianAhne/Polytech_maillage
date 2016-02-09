

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

typedef struct {
  int       np,na,nt,nr,ver,dim,mark;
  char     *namein,*nameout;
	double   *sol,o[3],rad;
  pPoint    point;
  pEdge     edge;
  pTria     tria;
} Mesh;
typedef Mesh * pMesh;


/* read mesh */
int loadMesh(pMesh mesh);
 
/** Save mesh data */
int saveMesh(pMesh mesh);

/* prototypes */
int saveSol(pMesh mesh,int it);



 
