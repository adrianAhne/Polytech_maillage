#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "libmesh5.h"
#include "mesh.h"
#include "rigidtransformation.h"
#include "distance.h"
#include "bucket.h"
#include "hash.h"
#include "ball.h" 
#include "hausdorff.h"
#include "InOut.h"



/* parse command line arguments (not used) */
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

 

int main(int argc,char *argv[]) {
 
	
	int test;

	fprintf(stdout,"*** Main3 (2016) ***\n ");

	fprintf(stdout,"TEST SELECTION\n\n Which test would you like to execute?\n\n 1. Rotation 2D\n 2. Rotation 3D\n 3. Superposition\n 4. Translation 2D\n 5. Translation 3D\n 6. Curvature 2D\n 7. Curvature 3D\n 8. Bucket \n 9. Distance point to triangle \n 10. Hash function \n 11. Normales \n 12. Distance point to mesh via bucket \n 13. Ball\n 14. distPointToBoule\n 15. distance Hausdorff\nChoose: \n");

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
		
		rotation2D(&mesh, angle);
		
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
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  		if ( !saveMesh(&mesh) )  return(1);
  		fprintf(stdout,"  -- WRITING COMPLETED\n The new mesh created is a (name).o.mesh \n ");	
	}
	
	
	/* Curvature 2D */
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
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  		if ( !saveSol(&mesh,1) )  return(1);
  		fprintf(stdout,"  -- WRITING COMPLETED\n \n ");		
	}
	
	
	/* Curvature 3D */
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
		
		/* save mesh */
		fprintf(stdout,"\n  -- OUTPUT DATA\n");
  		if ( !saveSol(&mesh,1) )  return(1);
  		fprintf(stdout,"  -- WRITING COMPLETED\n \n ");		
	}	


	/* Bucket */
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
		point.c[0] =   0.59149757132316072 ;
		point.c[1] =   0.52414378752244195 ;
		point.c[2] =   0.50863001183634415 ;
	
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
	   	int ret = bucket_retour_key( &bucket , &mesh, &point , 0.0 );
	   	fprintf(stdout," Retour de la fonction bucket = %d \n",ret); 
		free_bucket (&bucket);
				
	    fprintf(stdout,"  -- WRITING COMPLETED\n \n ");
	}

	
	/* Distance point to Triangle */
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


	/* Hash function */
	if (test == 10)
	{
		Mesh	mesh;
		int N,i;
		double result;

		double cb[3] = {.0,.0,.0};

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
		for (i = 0; i < 3; ++i)
			printf("Adj %d\n", mesh.adja[3*(5727-1)+1+i]/3);
			
		// Test de la fonction qui localise un point dans un triangle en brute force
		printf("Localise brute force: %d\n", localiseTriangleBruteForce(&mesh, &mesh.point[mesh.tria[5465].v[1]]));

		// Test de la fonction qui localise un point dans une triangle grace aux fonctions d'approche
		printf("locelt : %d\n", locelt(&mesh, 4493, &mesh.point[mesh.tria[5465].v[1]], cb));

		// Tester avec un point qui est situé au centre d'un triangle, calculer ses coordonnées cartésiennes
	}
	

	//Normales
 	if (test == 11)
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

		normalesOfTriangles(&mesh);
		mesh.nameout = "normales.mesh";
		saveMesh(&mesh);
 	}


	// distance using bucket
	if(test == 12 )
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

		int* VertToTria = (int*)calloc(mesh.np+1, sizeof(int));
		hashTria(&mesh, VertToTria);
		//printf("tab[2010]=%d\n", VertToTria[2010]);
		//printf("tab[2085]=%d\n", VertToTria[2085]);
		//printf("tab[2011]=%d\n", VertToTria[2011]);

		Point p;
		p.c[0] = .90;
		p.c[1] = .90;
		p.c[2] = .30;
		Bucket bucket ;
		int N;
		init_bucket( &bucket , &mesh ); 
		fill_bucket( &bucket , &mesh );
		fprintf(stdout,"\n  -- Creation bucket MESH \n\nPlease type the number of subdivision : \n");
		fflush(stdin);	// just to be able to do a scanf without conflicts
		fscanf(stdin,"%d",&N);
		bucket.size = N ;
		double dist;
		dist = distanceUsingBucket(&mesh, &p, VertToTria , &bucket);
		printf("DistanceUsingBucket = %f\n", dist);
		free(VertToTria);
 	}


 	// Bucket sort algorithm for the computation of the distance point-triangulation
 	if (test == 13)
 	{
 		Mesh	mesh;
		time_t debut;
		time_t fin;
		double difference;
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");

		
		/* default values */
		memset(&mesh,0,sizeof(Mesh));  	
		if ( !parsar(argc,argv,&mesh) )  return(1);

  		 
	  	/* read data */
  		fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		debut = time(NULL) ;
		Hedge *tab = (Hedge*)calloc(3*mesh.nt+1,sizeof(Hedge));
		hashHedge(&mesh, tab);
		setAdj(&mesh, tab);
		fin = time(NULL);
		difference = difftime (fin, debut);
		fprintf(stdout," temps d'éxécution %f \n",difference );
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");
		int** list = (int**)malloc(sizeof(int*)) ;
		int nb_triangle = boulep(&mesh,5727,0,list),i;
		printf("le point considéré = %d \n ", mesh.tria[5727].v[0] );
		
		
		printf("Nombre de triangles autour = %d \n ", nb_triangle );
		printf("Liste des triangles autour :  \n " );
		for (i=0;i<nb_triangle;i++)
			printf("triangle %d = %d \n ", i,(*list)[i] );
			
		printf("En utilisant les relations d'adjacences \n " );
		int** list2 = (int**)malloc(sizeof(int*)) ;
		int compt = boule_adj(&mesh, 5727, 0 , list2);
		printf("Nombre de triangles autour = %d \n ", compt );
		printf("Liste des triangles autour :  \n " );
		for (i=0;i<compt;i++)
			printf("triangle %d = %d \n ", i,(*list2)[i] );
		
		fprintf(stdout,"  -- WRITING COMPLETED\n \n ");
	}
		

	// test function hashTria
	if (test == 14)
	{
 		Mesh	mesh;
 		int i;
 		double dist;
		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");


		/* default values */
		memset(&mesh,0,sizeof(Mesh));  	
		if ( !parsar(argc,argv,&mesh) )  return(1);
  		 
  		/* read data */
  		fprintf(stdout,"\n  -- INPUT DATA MESH \n");
		if ( !loadMesh(&mesh) )  return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");	
		
	
		int* tab = (int*)calloc(mesh.np+1, sizeof(int));
		hashTria(&mesh, tab);
		printf("tab[877]=%d\n", tab[877]);
		//for(i=0; i<mesh.np+1; i++)
			//printf("i=%d triangle=%d\n", i, tab[i]);
		/*
		pTria tria = &mesh.tria[877];
		dist = distancePointToBoule(&mesh, 877, mesh.point[1377], tab);
		printf("dist=%f\n", dist);
		*/
			
		free(tab);
	}
 
 
 	// test hausdorff distance between two mesh
	if (test == 15)
	{
 		Mesh	meshA, meshB;

		/* parse arguments */
		fprintf(stdout,"\n  -- DATA MESH\n");

		/* default values */
		memset(&meshA,0,sizeof(Mesh));  	
		if ( !parsar(argc,argv,&meshA) )  return(1);
		if ( !loadMesh(&meshA) )  return(1);
		memset(&meshB,0,sizeof(Mesh));  	
		if ( !parsar(argc,argv,&meshB) )  return(1);
		if ( !loadMesh(&meshB) ) return(1);
		fprintf(stdout,"  -- DATA READING COMPLETED.\n");	
		
		double distHausdorff = Hausdorff(&meshA, &meshB);
		printf("distHausdorff = %f\n", distHausdorff);
 	}

	return(0);
}


