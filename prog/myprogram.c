#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libmesh5.h"
#include "mesh.h"
#include "groupfunctions.h"

#include "distanceMeshFunctions.h"

#include "bucket.h"
#include "hash.h"




/* read mesh */
int loadMesh(pMesh mesh) {
  pPoint     ppt;
  pEdge      pa;
  pTria      pt;
	double     xmin,ymin,zmin,xmax,ymax,zmax;
  float      fp1,fp2,fp3;
  int        k,i,inm,ia;
  char      *ptr,data[256];

  strcpy(data,mesh->namein);
	
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
	fprintf(stdout, " the xmax is %f \n ", xmax ) ;
	/* Center of the mesh*/
	mesh->o[0] = 0.5 * (xmin+xmax); 
	mesh->o[1] = 0.5 * (ymin+ymax); 
	mesh->o[2] = zmin; //0.5 * (zmin+zmax);
		
	//printf(" max %f %f  %f %f   %f %f\n",xmin,xmax,ymin,ymax,zmin,zmax);
	mesh->rad = max(xmax-xmin, max(ymax-ymin,zmax-zmin));
	mesh->rad *= 1.1;
	//printf("radius %f\n",mesh->rad);

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
  pNormal      pn;
  int          k,nr,outm;
  char         data[128];

	
  mesh->ver = GmfDouble;
  strcpy(data,mesh->nameout);
 fprintf(stdout,"ICI\n");
 fprintf(stdout,"Version = %d ; data = %s ; dimension = %d ", mesh->ver , data, mesh->dim );
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
	
// write normals 
if ( mesh->nn) {
  GmfSetKwd(outm,GmfNormals,mesh->nn);
  for(k=1; k<=mesh->nn; k++) {
    pn = &mesh->Normal[k];
    GmfSetLin(outm,GmfNormals,pn->n[0],pn->n[1],pn->n[2]);
  }
  GmfSetKwd(outm,GmfNormalAtVertices,mesh->np);
  for (k=1; k<=mesh->np; k++) {
    GmfSetLin(outm,GmfNormalAtVertices,k,k);
  }
}


  GmfCloseMesh(outm);
  return(1);
}

/* save solution */
int saveSol(pMesh mesh,int it) {
  double       buf[1];
  int          k,inm,type,typtab[GmfMaxTyp],is;
  char        *ptr,data[128];

	strcpy(data,mesh->nameout);
	// .o.mesh
	
	fprintf(stdout, " data = %s \n" , mesh->nameout ) ;  
  sprintf(data,"%s.sol",mesh->namein);

  if ( !(inm = GmfOpenMesh(data,GmfWrite,mesh->ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  mesh->ver = GmfDouble;
  type = 1;
  typtab[0] = GmfSca;

  GmfSetKwd(inm,GmfSolAtVertices,mesh->np,type,typtab);
  for (k=1; k<=mesh->np; k++) {     
    is = 3*(k-1)+1 ;
		buf[0] = mesh->sol[k];
    /*buf[1] = mesh->sol[is+1];
    buf[2] = mesh->sol[is+2];*/
    GmfSetLin(inm,GmfSolAtVertices,buf);
    //fprintf(stdout," solu = %f \n" , buf[0] ) ;
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



 

int main(int argc,char *argv[]) {
 
	
	int test;

  fprintf(stdout,"  -- Main3 (2016)\n");
	fprintf(stdout,"  TEST SELECTION\n\n Do you want to test ?: \n\n 1. Rotation 2D\n 2. Rotation 3D\n 3. Superposition\n 4. Translation 2D\n 5. Translation 3D\n 6. Courbure 2D\n 7. Courbure 3D\n 8. Bucket \n 9. Distance point to triangle \n 10. Hash function \n");
	fflush(stdin);
  fscanf(stdin,"%d",&test);
	
	
	/* Rotation 2D */
	if( test == 1 )
	{ 
		Mesh	mesh;
		float angle ;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Rotation 2D MESH \n\n Please type the angle for the rotation\n");
		fflush(stdin);
    fscanf(stdin,"%f",&angle);
		
		rotation2D(&mesh, angle)  ;
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveMesh(&mesh) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n The new mesh created is a (name).o.mesh \n");	
		
	}
	/* Rotation 3D */
	if(test==2)
	{
		Mesh	mesh;
		float angleX,angleY,angleZ ;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Rotation 2D MESH \n\n Please type the angles for the rotation\n angleX = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&angleX);
		
		fprintf(stdout,"\n  -- angleY = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&angleY);
		
		fprintf(stdout,"\n  -- angleZ = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&angleZ);
		
		
		rotation3D(&mesh, angleX, angleY, angleZ);
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveMesh(&mesh) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n The new mesh created is a (name).o.mesh \n ");	
	}
	/* Superposition */
	if( test == 3 )
	{
		Mesh	mesh1,mesh2,mesh3;
		
		/* default values */
		memset(&mesh1,0,sizeof(Mesh));
		memset(&mesh2,0,sizeof(Mesh));
		memset(&mesh3,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH1\n");
  	if ( !parsar(argc,argv,&mesh1) )  return(1);
  	fprintf(stdout,"\n  -- DATA MESH2\n");
  	if ( !parsar(argc,argv,&mesh2) )  return(1);
  	fprintf(stdout,"\n  -- NAME FOR THE NEW MESH \n");
  	if ( !parsar(argc,argv,&mesh3) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH1 \n");
		if ( !loadMesh(&mesh1) )  return(1);
		fprintf(stdout,"\n  -- INPUT DATA MESH2 \n");
		if ( !loadMesh(&mesh2) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Superposition MESH \n");
		Superposition(&mesh1,&mesh2 , &mesh3 ) ;
		if ( ! (&mesh3)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveMesh(&mesh3) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n ");	
	}
	
	/* Translation 2D */
	if( test == 4 )
	{
		Mesh	mesh;
		float lengthX,lengthY;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Translation 2D MESH \n\n Please type the lengths for the translation\n lengthX = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&lengthX);
		
		fprintf(stdout,"\n  -- lengthY = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&lengthY);

		translation2D(&mesh, lengthX,lengthY);
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveMesh(&mesh) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n The new mesh created is a (name).o.mesh \n ");	
	} 
	/* Translation3D */ 
	if (test == 5 )
	{
		Mesh	mesh;
		float lengthX,lengthY,lengthZ;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Translation 3D MESH \n\n Please type the lengths for the translation\n lengthX = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&lengthX);
		
		fprintf(stdout,"\n  -- lengthY = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&lengthY);
    
    fprintf(stdout,"\n  -- lengthZ = ? \n");
		fflush(stdin);
    fscanf(stdin,"%f",&lengthZ);

		translation3D(&mesh, lengthX,lengthY,lengthZ);
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveMesh(&mesh) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n The new mesh created is a (name).o.mesh \n ");	
	}
	/* Courbure 2D */
	if ( test == 6 )
	{
		Mesh	mesh;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Courbure 2D MESH \n\n ");

		if ( ! courbure2D(&mesh)) return(1);
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveSol(&mesh,1) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n \n ");		
	}
	/* Courbure 3D */
	if ( test == 7 )
	{
		Mesh	mesh;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		/* FUNCTION */
		fprintf(stdout,"\n  -- Courbure 3D MESH \n\n ");

		if ( ! courbure3D(&mesh)) return(1);
		if ( ! (&mesh)) return(1);
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  	if ( !saveSol(&mesh,1) )  return(1);
  	fprintf(stdout,"  -- WRITING COMPLETED\n \n ");		
	}	
if ( test == 8 ) 
	{
		Mesh	mesh;
		Bucket bucket;
		int N,indice =5207;
    Point point;
		               
	
		
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  	if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  	 /* read data */
  	fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		
		                                        
			/* init point */            
		point.c[0] =  -3.29169692485065  ;
		point.c[1] =  -1.6280599609887  ;
		point.c[2] =  1.58560780723174  ;
	
		
		fprintf(stdout," pointx = %f \n " , mesh.point[5207].c[0] ) ;
    fprintf(stdout," pointx = %f \n " , mesh.point[5207].c[1] ) ;
    fprintf(stdout," pointx = %f \n " , mesh.point[5207].c[2] ) ;
		
		positive_boundingbox( &mesh , &point );
		/* FUNCTION */
		fprintf(stdout,"\n  -- Creation bucket MESH \n\n Please type the number of subdivision : \n");
		fflush(stdin);
    fscanf(stdin,"%d",&N);
    bucket.size = N ;
    printf(" size = %d ", bucket.size);
    int i = max(0,(int)(N*point.c[0])-1) ;
		int j = max(0,(int)(N*point.c[1])-1) ;
		int k = max(0,(int)(N*point.c[2])-1) ;
    fprintf(stdout," i = %d p = %f  \n " , i, point.c[0]) ;
		fprintf(stdout," j = %d p = %f  \n " , j , point.c[1]) ;
		fprintf(stdout," k = %d  p = %f \n " , k , point.c[2]) ;
    fprintf(stdout, " key point = %d \n", (j*N+k)*N+i) ;
    
		init_bucket( &bucket , &mesh); 
		fill_bucket( &bucket , &mesh ) ;
    use_bucket( &bucket , &mesh, &point , 0.0 );
		free_bucket (&bucket);
		
		//if ( ! (&mesh)) return(1);
    
		
    fprintf(stdout,"  -- WRITING COMPLETED\n \n ");
	}

	
	if ( test == 9 )
	{
		Mesh	mesh;
		int N;
		double result;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  		if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  		 /* read data */
  		fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");

		result = distPointToTriangle(&mesh, &mesh.tria[1], &mesh.point[mesh.tria[1].v[0]]);
		printf("%f\n", result);
	}

	if (test == 10)
	{
		Mesh	mesh;
		int N,i;
		double result;
		/* default values */
		memset(&mesh,0,sizeof(Mesh));
		
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");
  		if ( !parsar(argc,argv,&mesh) )  return(1);
  	
  		 /* read data */
  		fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");

		Hedge *tab = (Hedge*)calloc(3*mesh.nt+1,sizeof(Hedge));
		hashHedge(&mesh, tab);
		setAdj(&mesh, tab);
		// Test de la fonction qui localise un point dans un triangle en brute force
		printf("%d\n", localiseTriangleBruteForce(&mesh, &mesh.point[mesh.tria[15].v[1]]));
		

	}

 

 

	


	
	/*if ( !courbure2D(&mesh) ) return (1) ;
	if ( !saveSol(&mesh,1) )  return(1);
	*/



	/*rotation3D(&mesh, 3.14 , 0 , 0);*/

	//Adrian: Translation 2D
	//translation2D(&mesh, 0.5, 0);
	//translation3D(&mesh, 0, 1, 1);
		

  // Tupac : write data into a new file




  return(0);
}
