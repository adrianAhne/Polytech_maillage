#ifndef HASH_H
#define HASH_H

/*
Hedge struct doc:
na et nb les extrémités de l'arête

ia = min(na,nb)
ib = max(na,nb)

Ces objects Hedge sont rangés dans un tableau global

key = (KA*min(na,nb)+KB*max(na,nb))%hsize
KA = 7, KB = 11, hsize = mesh->np
*/

typedef struct {
	int ia, ib, adj1, adj2, nxt;
} Hedge;
typedef Hedge * pHedge;

int hashHedge(pMesh mesh);

int setAdj(pMesh mesh);


#endif