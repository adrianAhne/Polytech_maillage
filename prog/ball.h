#ifndef BALL_H
#define BALL_H



/* Searches the ball/triangles around the given point and saves them in list
   Brute force method, walking trough all triangles and check if points is 
   attached or not
   triangles are saved like this: list[i] = 3*k+j with k is the k-th triangle of the ball and j the position of the point in the ball (0,1,2) */
int boulep(pMesh mesh, int start, int point , int** list);

/* Function boule optimised!
   Searches the ball/triangles around the given point and saves them in list
   Does not walk anymore through all the triangles. Checks just in the 
   neighbourhood if the triangles belong to the point
   triangles are saved like this: list[i] = 3*k+j with k is the k-th triangle of the ball and j the position of the point in the ball (0,1,2) */
int boule_adj(pMesh mesh, int start, int point , int** list);

/* This function takes a array tab of the size of all points in the mesh
 For each point in the mesh it associates one possible triangle belonging to this point */
void hashTria(pMesh mesh, int *tab);

#endif

