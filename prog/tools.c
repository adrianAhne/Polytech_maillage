#include "mesh.h"

extern Info  info;




/* Compute Hausdorff distance between two set of triangles discarding bounding box and interior Dirichlet boundary */
/* mesh1 : target mesh, mesh2 : template mesh */
double hausdorff_bruteforce(pMesh mesh1, pMesh mesh2){
  pTria   pt,pt1;
  pPoint  p0,p1,pmil,p2,p3,p4,p5;
  Point   mil;
  double  rho2,haus,d, d0,d1,d2,dmil;
  double rho1;
  int     k,i,j,nac1,nac2;
  char    proj;
  
  rho1 = 0.0;
  rho2 = 0.0;
  pmil = &mil;
  
  nac1 = 0;
  nac2 = 0;
  
  for(k=1; k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    if ( pt->ref != mesh2->refb)
    pt->ref = 0;
  }
  
  for(k=1; k<=mesh2->nt;k++){
    pt1 = &mesh2->tria[k];
    if ( pt1->ref != mesh2->refb)
    pt1->ref = 0;
  }
  
  /* Set active triangles for mesh 1 (discard bounding box) */
  for(k = 1; k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    p0 = &mesh1->point[pt->v[0]];
    p1 = &mesh1->point[pt->v[1]];
    p2 = &mesh1->point[pt->v[2]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)||(p0->c[1]<0.01)||(p0->c[1]>0.99)||(p0->c[2]<0.01)||(p0->c[2]>0.99)) continue;
    pt->ref  =5;
    nac1++;
  }
  
  /* Set active triangles for mesh 2 (discarding Dirichlet) */
  for(k = 1; k<=mesh2->nt;k++){
    pt = &mesh2->tria[k];
    if ( pt->ref != mesh2->refb) {
      pt->ref  =5;
      nac2++;
    }
  }

  printf("Number of active triangles %d %d \n", nac1, nac2);

  /* Compute rho(\Gamma_1,\Gamma_2)*/
  for(k=1;k<=mesh1->nt;k++){
    pt = &mesh1->tria[k];
    if(pt->ref ==5){
    
    p0 = &mesh1->point[pt->v[0]];
    p1 = &mesh1->point[pt->v[1]];
    p2 = &mesh1->point[pt->v[2]];
      
    for (i=0; i<3;i++) pmil->c[i] = (p0->c[i]+p1->c[i]+p2->c[i])/3;

    d0 = 10.0;
    d1 = 10.0;
    d2 = 10.0;
    dmil = 10.0;
      
    
    for(j=1;j<=mesh2->nt;j++){
      pt1 = &mesh2->tria[j];
      if ( pt1->ref ==5 ) {
      p3 = &mesh2->point[pt1->v[0]];
      p4 = &mesh2->point[pt1->v[1]];
      p5 = &mesh2->point[pt1->v[2]];
      
      d = distpt_3d(p3,p4,p5,p0,&proj);
      d0 = LS_MIN(d0,d);
      
      d = distpt_3d(p3,p4,p5,p1,&proj);
      d1 = LS_MIN(d1,d);
      
      d = distpt_3d(p3,p4,p5,p2,&proj);
      d2 = LS_MIN(d2,d);
      
      d = distpt_3d(p3,p4,p5,pmil,&proj);
      dmil = LS_MIN(dmil,d);
      
         }
       }
    }
    rho1 = LS_MAX(rho1,d0);
    rho1 = LS_MAX(rho1,d1);
    rho1 = LS_MAX(rho1,d2);
  }
  printf("rho1 %e  \n", rho1);
    
    /* Compute rho(\Gamma_1,\Gamma_2)*/
    for(k=1;k<=mesh2->nt;k++){
      pt = &mesh2->tria[k];
      if(pt->ref ==5){
        
        p0 = &mesh2->point[pt->v[0]];
        p1 = &mesh2->point[pt->v[1]];
        p2 = &mesh2->point[pt->v[2]];
        
        for (i=0; i<3;i++) pmil->c[i] = (p0->c[i]+p1->c[i]+p2->c[i])/3;
        
        d0 = 10.0;
        d1 = 10.0;
        d2 = 10.0;
        dmil = 10.0;
        
        
        for(j=1;j<=mesh1->nt;j++){
          pt1 = &mesh1->tria[j];
          p3 = &mesh1->point[pt1->v[0]];
          p4 = &mesh1->point[pt1->v[1]];
          p5 = &mesh1->point[pt1->v[2]];
          
          d = distpt_3d(p3,p4,p5,p0,&proj);
          d0 = LS_MIN(d0,d);
          
          d = distpt_3d(p3,p4,p5,p1,&proj);
          d1 = LS_MIN(d1,d);
          
          d = distpt_3d(p3,p4,p5,p2,&proj);
          d2 = LS_MIN(d2,d);
          
          d = distpt_3d(p3,p4,p5,pmil,&proj);
          dmil = LS_MIN(dmil,d);
          
          
        }
      }
      rho2 = LS_MAX(rho1,d0);
      rho2 = LS_MAX(rho1,d1);
      rho2 = LS_MAX(rho1,d2);
      rho2 = LS_MAX(rho1,dmil);
    }
  
  haus = LS_MAX(rho1,rho2);
  haus = sqrt(haus);
  
  return(haus);
}
  

  /* Compute squared distance from pq to tria p0p1p2,
   proj = 1: projection onto face, =2: distance to vertex or edge */
double distpt_3d(pPoint p0,pPoint p1,pPoint p2,pPoint pq,char *proj) {
    Point pointPlan0, pointPlan1, pointPlan2, pointPlana;
    double lx1, ly1, lz1, lx2, ly2, lz2, lxq, lyq, lzq, longp0p1, cosAlpha, sinAlpha, aa, bb, ab, ll, l;
    double p0XTemp, p0YTemp, p0ZTemp, p1XTemp, p1YTemp, p1ZTemp, p2XTemp, p2YTemp, p2ZTemp, pqXTemp, pqYTemp, pqZTemp; // pour stocker temporairement
    double m11, m12, m13, m21, m22, m23, m31, m32, m33, d01, d12, d02, dTmp;
    double zone01, zone02, zone12;
    
    *proj = 1;
    lx1 = p1->c[0] - p0->c[0];
    ly1 = p1->c[1] - p0->c[1];
    lz1 = p1->c[2] - p0->c[2];
    
    lx2 = p2->c[0] - p0->c[0];
    ly2 = p2->c[1] - p0->c[1];
    lz2 = p2->c[2] - p0->c[2];
    
    lxq = pq->c[0] - p0->c[0];
    lyq = pq->c[1] - p0->c[1];
    lzq = pq->c[2] - p0->c[2];
    
    longp0p1 = sqrt(lx1*lx1 + ly1*ly1 + lz1*lz1);
    
    cosAlpha = lz1/longp0p1;
    sinAlpha = sqrt(1.0-cosAlpha*cosAlpha);
    
    p0XTemp = 0.0;
    p0YTemp = 0.0;
    p0ZTemp = 0.0;
    p1XTemp = lx1;
    p1YTemp = ly1;
    p1ZTemp = lz1;
    p2XTemp = lx2;
    p2YTemp = ly2;
    p2ZTemp = lz2;
    pqXTemp = lxq;
    pqYTemp = lyq;
    pqZTemp = lzq;
    
    /* Apply rotation with angle alpha */
    
    aa = lx1*lx1;
    bb = ly1*ly1;
    ab = lx1*ly1;
    ll = aa + bb;
    l = sqrt(ll);
    
    if(ll<EPS1){
      if(p1ZTemp <= 0.0){
        p1ZTemp = -p1ZTemp;
        p2ZTemp = -p2ZTemp;
        pqZTemp = -pqZTemp;
      }
    }
    
    else{
      m11 = (aa*cosAlpha + bb)/ll;
      m12 = (ab*cosAlpha - ab )/ll;
      m13 = -lx1*sinAlpha/l;
      
      m21 = (ab * cosAlpha - ab)/ll;
      m22 = (bb*cosAlpha + aa)/ll;
      m23 = - ly1 * sinAlpha/l;
      
      m31 = (lx1* sinAlpha)/l;
      m32 = (ly1 * sinAlpha)/l;
      m33 = cosAlpha;
      
      p1XTemp = m11*lx1 + m12*ly1 + m13 * lz1;
      p1YTemp = m21*lx1 + m22*ly1 + m23 * lz1;
      p1ZTemp = longp0p1;
      
      p2XTemp = m11*lx2 + m12*ly2 + m13 * lz2;
      p2YTemp = m21*lx2 + m22*ly2 + m23 * lz2;
      p2ZTemp = m31*lx2 + m32*ly2 + m33 * lz2;
      
      pqXTemp = m11*lxq + m12*lyq + m13 * lzq;
      pqYTemp = m21*lxq + m22*lyq + m23 * lzq;
      pqZTemp = m31*lxq + m32*lyq + m33 * lzq;
    }
    
    /* Apply second rotation to put point p2 in plane (yz), in half plane y > 0*/
    assert((p2XTemp* p2XTemp + p2YTemp* p2YTemp) > 0.0);  // au cas ou...
    cosAlpha = p2YTemp/(sqrt(p2XTemp* p2XTemp + p2YTemp* p2YTemp));
    sinAlpha = sqrt(1.0-cosAlpha * cosAlpha);
    if(p2XTemp <=0.0) {sinAlpha = -sinAlpha; }
    
    lx2 = 0.0;
    ly2 = sinAlpha * p2XTemp + cosAlpha * p2YTemp;
    lz2 = p2ZTemp;
    
    lxq = cosAlpha * pqXTemp - sinAlpha * pqYTemp;
    lyq = sinAlpha * pqXTemp + cosAlpha * pqYTemp;
    lzq = pqZTemp;
    
    zone01 = lyq;
    zone02 = ly2*lzq - lz2*lyq;
    zone12 = lyq*(lz2-longp0p1) - ly2*(lzq - longp0p1);
    
    if((zone01>=0.0)&&(zone02 >=0.0)&&(zone12>=0.0))
    {
      return (lxq*lxq);
    }
    
    else
    {
      pointPlan0.c[0] = 0.0;
      pointPlan0.c[1] = 0.0;
      
      pointPlan1.c[0] = 0.0;
      pointPlan1.c[1] = longp0p1;
      
      pointPlan2.c[0] = ly2;
      pointPlan2.c[1] = lz2;
      
      pointPlana.c[0] = lyq;
      pointPlana.c[1] = lzq;
      
      d01 = distpt_23d(&pointPlan0, &pointPlan1, &pointPlana);
      d02 = distpt_23d(&pointPlan0, &pointPlan2, &pointPlana);
      d12 = distpt_23d(&pointPlan1, &pointPlan2, &pointPlana);
      
      dTmp = LS_MIN(d01, LS_MIN(d02, d12));
      *proj = 2;
      
      return(dTmp + lxq*lxq);
    }
    
  }
