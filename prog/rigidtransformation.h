#ifndef RIGIDTRANSFORMATION_H
#define RIGIDTRANSFORMATION_H

#define PI         3.14159265

// translates every point in the mesh by lengthX in x-direction and lengthY in ydirection in the 2 dimensional case
void translation2D(pMesh mesh, double lengthX, double lengthY);

// translates every point in the mesh by lengthX in x-direction and lengthY in y-direction and lengthZ in z-direction in the 3 dimensional case
void translation3D(pMesh mesh, double lengthX, double lengthY, double lengthZ);


/* rotation via the angle angle in the two dimension case */
void rotation2D(pMesh mesh, double angle);

/* rotation via the angles angleX, angleY, angleZ in the three dimension case */
void rotation3D(pMesh mesh, double angleX, double angleY, double angleZ);


/* Changes from 2D to 3D by setting the z-coordinate to 0 */
void Change2Dto3D( pMesh Mesh );

/* Combines all information (points, edges, triangles,...) of two meshes in a new mesh*/
int Superposition(pMesh Mesh1, pMesh Mesh2, pMesh Mesh_final ) ;



#endif
