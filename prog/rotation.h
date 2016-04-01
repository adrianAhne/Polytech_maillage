#ifndef ROTATION_H
#define ROTATION_H

#define PI         3.14159265

/* FUNCTION rotation2D 
This function will rotate a 2D mesh for a given angle. 
PARAMETERS : a pointer through the mesh and the angle
*/
void rotation2D(pMesh mesh, double angle);

/*FUNCTION rotation3D
This function will rotate a 3D mesh. You can give angles for x,y,z axis.
PARAMETERS : a pointer through a mesh, 3 angles */
void rotation3D(pMesh mesh, double angleX, double angleY, double angleZ);

#endif

