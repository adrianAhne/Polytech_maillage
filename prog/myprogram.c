#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"




/* read mesh */
int loadMesh(pMesh mesh) {
  pPoint     ppt;
  pEdge      pa;
  pTria      pt;
	double     xmin,ymin,zmin,xmax,ymax,zmax;
  float      fp1,fp2,fp3;
  int        k,i,inm,ia;
  char      *ptr,data[256];
	fprintf(stdout,"ici\n");
  strcpy(data,mesh->namein);
	fprintf(stdout,"là\n");
	/* test if we have a ".mesh"*/
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->na = GmfStatKwd(inm,GmfEdges);
  mesh->nr = GmfStatKwd(inm,GmfRidges);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);

  if ( !mesh->np || !mesh->nt ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }




  /* memory allocation */
  mesh->point = (pPoint)calloc(mesh->np+1,sizeof(Point));
  assert(mesh->point);
	mesh->sol = (double*)calloc(3*mesh->np+1, sizeof(double));
	assert(mesh->sol);
  if ( mesh->na ){
    mesh->edge = (pEdge)calloc(mesh->na+1,sizeof(Edge));
    assert(mesh->edge);
  }
  if ( mesh->nt ) {
    mesh->tria = (pTria)calloc(mesh->nt+1,sizeof(Tria));
    assert(mesh->tria);
  }


	/* fill the tab of points */

  GmfGotoKwd(inm,GmfVertices);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( mesh->ver == GmfFloat ) {
      GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
      ppt->c[0] = fp1;
      ppt->c[1] = fp2;
      ppt->c[2] = fp3;
    }
    else
      GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
		
  }

	/* Center of the mesh*/
	mesh->o[0] = 0.5 * (xmin+xmax); 
	mesh->o[1] = 0.5 * (ymin+ymax); 
	mesh->o[2] = zmin; //0.5 * (zmin+zmax);
		
	printf(" max %f %f  %f %f   %f %f\n",xmin,xmax,ymin,ymax,zmin,zmax);
	mesh->rad = max(xmax-xmin, max(ymax-ymin,zmax-zmin));
	mesh->rad *= 1.1;
	printf("radius %f\n",mesh->rad);

  /* read triangles and fill the tab */
  GmfGotoKwd(inm,GmfTriangles);
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    GmfGetLin(inm,GmfTriangles,&pt->v[0],&pt->v[1],&pt->v[2],&pt->ref);
		for (i=0; i<3; i++)
			mesh->point[pt->v[i]].ref = pt->ref;
  }

  /* read edges */
  GmfGotoKwd(inm,GmfEdges);
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
    GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
  }
	if ( mesh->nr ) {
	  GmfGotoKwd(inm,GmfRidges);
    for (k=1; k<=mesh->nr; k++) {
      GmfGetLin(inm,GmfRidges,&ia);
			pa = &mesh->edge[ia];
			pa->tag = 1;
	  }
  }

  fprintf(stdout,"  %%%% NUMBER OF VERTICES   %8d\n",mesh->np);
  fprintf(stdout,"  %%%% NUMBER OF TRIANGLES  %8d\n",mesh->nt);  
  fprintf(stdout,"  %%%% NUMBER OF EDGES      %8d (%d ridges)\n",mesh->na,mesh->nr);
	fprintf(stdout,"  %%%% DIMENSION  %8d\n",mesh->dim);  
	fprintf(stdout,"  %%%% MARK  %8d\n",mesh->mark);  
	fprintf(stdout,"  %%%% VER  %8d\n",mesh->ver); 

  GmfCloseMesh(inm);
  return(1);
}


/** Save mesh data */
int saveMesh(pMesh mesh) {
  pPoint       ppt;
  pTria        pt;
  pEdge        pa;
  int          k,nr,outm;
  char         data[128];

	
  mesh->ver = GmfDouble;
  strcpy(data,mesh->nameout);
 fprintf(stdout,"ICI\n");
  if ( !(outm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
 
  /* vertices */
  GmfSetKwd(outm,GmfVertices,mesh->np);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    GmfSetLin(outm,GmfVertices,ppt->c[0],ppt->c[1],ppt->c[2],ppt->ref);
  }

  GmfSetKwd(outm,GmfTriangles,mesh->nt);
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    GmfSetLin(outm,GmfTriangles,pt->v[0],pt->v[1],pt->v[2],pt->ref);
  }

  if ( mesh->na ) {
    GmfSetKwd(outm,GmfEdges,mesh->na);
    for (k=0; k<=mesh->na; k++) {
			pa = &mesh->edge[k];
			GmfSetLin(outm,GmfEdges,pa->v[0],pa->v[1],pa->ref);
    }
    if ( mesh->nr ) {
      GmfSetKwd(outm,GmfRidges,mesh->nr);
      nr = 0;
      for (k=0; k<=mesh->na; k++) {
        pa = &mesh->edge[k];
        if ( !pa->tag )  continue;
				nr++;
        GmfSetLin(outm,GmfRidges,nr);
      }
    }
	}

  GmfCloseMesh(outm);
  return(1);
}

/* save solution */
int saveSol(pMesh mesh,int it) {
  double       buf[3];
  int          k,inm,type,typtab[GmfMaxTyp],is;
  char        *ptr,data[128];

  strcpy(data,mesh->nameout);
  sprintf(data,"%s.sol",mesh->nameout);

  if ( !(inm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  mesh->ver = GmfDouble;
  type = 1;
  typtab[0] = GmfVec;

  GmfSetKwd(inm,GmfSolAtVertices,mesh->np,type,typtab);
  for (k=1; k<=mesh->np; k++) {     
    is = 3*(k-1)+1;
		buf[0] = mesh->sol[is+0];
		buf[1] = mesh->sol[is+1];
		buf[2] = mesh->sol[is+2];
    GmfSetLin(inm,GmfSolAtVertices,buf);
  }
  GmfCloseMesh(inm);
  return(1);
}


/* parse command line arguments */
static int parsar(int argc,char *argv[],pMesh mesh) {
  int     i;
  char   *ptr;

  i = 1;
  while ( i < argc ) {
    if ( mesh->namein == NULL ) {
      mesh->namein = (char*) calloc(strlen(argv[i])+1,sizeof(char));
      strcpy(mesh->namein,argv[i]);
    }
    else if ( mesh->nameout == NULL ){
      mesh->nameout = (char*) calloc(strlen(argv[i])+1,sizeof(char));
      strcpy(mesh->nameout,argv[i]);
    }
    i++;
  }

  /* check file names */
  if ( mesh->namein == NULL ) {
    mesh->namein = (char *)calloc(128,sizeof(char));
    assert(mesh->namein);
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh->namein);
  }
  if ( mesh->nameout == NULL ) {
    mesh->nameout = (char *)calloc(128,sizeof(char));
    assert(mesh->nameout);
    strcpy(mesh->nameout,mesh->namein);
    ptr = strstr(mesh->nameout,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(mesh->nameout,".o.mesh");
    ptr = strstr(mesh->nameout,".meshb");
    if ( ptr )  strcat(mesh->nameout,"b");
  }

  return(1);
}


int hello(pMesh mesh,int ref) {
  
  fprintf(stdout,"  -- HELLO ref = %d \n",ref);
  
  return(1);
}
/*** FUNCTION SUPERPOSITION ***/
/*	Parameters : 3 meshs: the 2 meshs to combine + the combinaison of the 2
	return : 1 if it works, 0 if not
	This function will join a mesh to another 
	Before this function The function loadMesh has to be done on the two meshs to combine, 
	in order to fill the number of Vertices, Triangles and edges of each one.
	After that we will go through the file of Mesh1 as a text file and find where is the last point
*/

int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) 
{
	
	/* Variables */
	int i,j;
	
	/* First we put the number of vertices, triangles and edges of the final mesh */ 
	Mesh_final->nt = Mesh1->nt + Mesh2->nt 	;
	Mesh_final->np = Mesh1->np + Mesh2->np 	;
	Mesh_final->na = Mesh1->na + Mesh2->na 	;
	Mesh_final->nr = Mesh1->nr + Mesh2->nr 	;

	/* Fill the dimension ( we take the more important dimension) and the mark */

	Mesh_final->dim = Mesh1->dim ;
	Mesh_final->mark = Mesh1->mark ;
	
	
	/* Version of the mesh wanted */
	Mesh_final->ver = 1 ;
	

	/* Memory allocation */

	Mesh_final->point = (pPoint)calloc(Mesh_final->np,sizeof(Point));
  assert(Mesh_final->point);
	Mesh_final->sol = (double*)calloc(3*Mesh_final->np, sizeof(double));
	assert(Mesh_final->sol);
  if ( Mesh_final->na ){
    Mesh_final->edge = (pEdge)calloc(Mesh_final->na,sizeof(Edge));
    assert(Mesh_final->edge);
  }
  if ( Mesh_final->nt ) {
    Mesh_final->tria = (pTria)calloc(Mesh_final->nt,sizeof(Tria));
    assert(Mesh_final->tria);
  }
	
	/* Now we fill the tab of points, vertices and edges */

	/* For the tab of points we first input in the new tab the vertices of mesh1. 
		 After that we eliminate the last char of the tab ("\0") and input the vertices of the mesh2
	*/  
	
	
		for(i=0;i<=Mesh1->np;i++)
		{
			Mesh_final->point[i] = Mesh1->point[i];
		}
	
	/*points of mesh2*/
	
		for(j=0;j<= Mesh2->np;j++)
		{
			Mesh_final->point[i+j] = Mesh2->point[j];
		}
	
	/* Triangles of mesh1 */
	
		for(i=0;i<=Mesh1->nt;i++)
		{
			Mesh_final->tria[i] = Mesh1->tria[i];
		}

	/*Triangles of mesh2*/
	
		for(j=0;j<= Mesh2->nt;j++)
		{
			Mesh_final->tria[i+j].v[0] = (Mesh2->tria[j+1].v[0])+((Mesh1->np) + 1 );
			Mesh_final->tria[i+j].v[1] = (Mesh2->tria[j+1].v[1])+((Mesh1->np) + 1);
			Mesh_final->tria[i+j].v[2] = (Mesh2->tria[j+1].v[2])+((Mesh1->np) + 1);
		}


	/* the center of the mesh is the average of the center of the 2 meshs */


	
	

	

	return(1);
}


/***	FUNCTION ROTATION_2D	***/
/*	Parameters: a mesh
	Return: the mesh rotated
	We will go through all the points of the mesh and for each one we will apply the rotation and write the new coordinates in a new mesh file
*/

 

int main(int argc,char *argv[]) {
  Mesh  mesh1;
	Mesh	mesh2;
	Mesh	mesh3;

  fprintf(stdout,"  -- Main3 (2016)\n");

  /* default values */
  memset(&mesh1,0,sizeof(Mesh));
	memset(&mesh2,0,sizeof(Mesh));
	memset(&mesh3,0,sizeof(Mesh));

  /* parse arguments */
	fprintf(stdout,"\n  -- DATA MESH1\n");
  if ( !parsar(argc,argv,&mesh1) )  return(1);

	fprintf(stdout,"\n  -- DATA MESH2\n");
  if ( !parsar(argc,argv,&mesh2) )  return(1);

	fprintf(stdout,"\n  -- DATA MESH3\n");
  if ( !parsar(argc,argv,&mesh3) )  return(1);
  /* read data */

  fprintf(stdout,"\n  -- INPUT DATA MESH1 \n");
  if ( !loadMesh(&mesh1) )  return(1);
  fprintf(stdout,"  -- DATA READING COMPLETED.\n");

	fprintf(stdout,"\n  -- INPUT DATA MESH2 \n");
  if ( !loadMesh(&mesh2) )  return(1);
  fprintf(stdout,"  -- DATA READING COMPLETED.\n");
	

	if ( ! (Superposition(&mesh1, &mesh2, &mesh3 ) )) return(1) ;

	if ( ! (&mesh3)) return(1);


	

  if ( !saveMesh(&mesh3))  return(1);
  if ( !hello(&mesh3,2) )  return(1);

  fprintf(stdout,"  -- WRITING COMPLETED\n");
  return(0);
}
