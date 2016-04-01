#ifndef TRANSLATION_H
#define TRANSLATION_H

/* FUNCTION translation 2D 
For this function we will take the lenght in X and Y directions and translate the mesh. 
In order to do that we go through all the point and change their coordinates.
PARAMETERS : a pointer through the mesh and the given lengths (x,y) for the modification 
*/
void translation2D(pMesh mesh, double lengthX, double lengthY);

/* FUNCTION translation 2D 
Same idea as the 2D translation 
PARAMETERS : a pointer through the mesh and the given lengths (x,y,z)  for the modification 
*/
void translation3D(pMesh mesh, double lengthX, double lengthY, double lengthZ);

#endif
