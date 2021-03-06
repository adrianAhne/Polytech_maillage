#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "libmesh5.h"
#include "mesh.h"



/* 
	Contains all the functions who are responsable to read the datas, to load them, to save them,...
	They use mainly the functions defined in the library libmesh5
*/




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
	
	/* test if we searched file .mesh exists */
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

	//  reading information out of file and saves them in the mesh
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

	// fprintf(stdout, " the xmax is %f \n ", xmax ) ;
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
    	// add to structure point an integer which indicates a triangle he belongs to
		for (i=0; i<3; i++){
			ppt = &mesh->point[pt->v[i]];	
			if(!(ppt->s))
				ppt->s = k;	
		}
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


/* Save mesh data */
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

