#ifndef INOUT_H
#define INOUT_H

/* read mesh */
int loadMesh(pMesh mesh);

/** Save mesh data */
int saveMesh(pMesh mesh);

/* save solution */
int saveSol(pMesh mesh,int it);

#endif
