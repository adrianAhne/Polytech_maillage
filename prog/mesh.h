

#define NPMAX 10000

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/* Structure point */
typedef struct {
  double    c[3];
  int ref, s;
} Point;
typedef Point * pPoint;

/* Structure normal */
typedef struct {
  double    n[3];
} Normal;
typedef Normal * pNormal; 

/* Structure edge */
typedef struct {
  int       v[2],ref,tag;
} Edge;
typedef Edge * pEdge;
  
/* Structure triangle */
typedef struct {
  int       v[3],ref;
} Tria;
typedef Tria * pTria;

/* Structure normales on the triangle  */
typedef struct {
  // Normal vertice defined by 2 points : norm start and norm end
  double  n[3];
  double  weight;
} TriaNorm;
typedef TriaNorm * pTriaNorm;

/* Main structure describing the mesh */
typedef struct {
  int       np,na,nt,nr,nn,ver,dim,mark; // number of points, number of triangles, dimension,...
  char     *namein,*nameout;
	double   *sol,o[3],rad;
  pPoint    point;
  pNormal   Normal;
  pEdge     edge;
  pTria     tria;
  pTriaNorm triaNorm;
  int *adja;	// adjacents
} Mesh;
typedef Mesh * pMesh;


/* read mesh */
//int loadMesh(pMesh mesh);
 
/** Save mesh data */
//int saveMesh(pMesh mesh);

/* prototypes */
//int saveSol(pMesh mesh,int it);



 

