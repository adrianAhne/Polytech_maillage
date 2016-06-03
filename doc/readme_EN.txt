******************* Start the program ****************************

Steps:

A. Go into the directory ..Polytech_maillage/prog

B. Compile the program with the command "make"

C. Start the program with the command "./main"

	Maintenant vous pouvez choisir entre les cas de test suivants:
	
	1. Rotation 2D
	   rotation via the angle 'angle' in dimension 2
	   Input:	- mesh mesh
	   			- angle
	   Output:	- points of the mesh are turned
	   
	2. Rotation 3D
	   rotation via the angles angleX, angleY, angleZ in dimension 3
	   Input:	- mesh mesh
	   			- angle
	   Output:	- points of the mesh are turned
	   			
	3. Superposition
	   Fusion of all the information (points, edge, triangles,...) of the two
	   mesh Mesh1 et Mesh2 in a new mesh Mesh_final
	   Input:	- mesh Mesh1, Mesh2, Mesh_final)
	   Output:  - all information are saved in the new mesh 'Mesh_final'
	   
	4. Translation 2D
	   translate each point of the mesh by the values lengthX, lengthY
	   in dimension 2
	   Input: 	- mesh mesh
	   			- lengthX, lengthY
	   Output:	- points of the mesh are translated
	   			
	5. Translation 3D
	   translate each point of the mesh by the values lengthX, lengthY,
	   lengthZ in dimension 3
	   Input: 	- mesh mesh
	   			- lengthX, lengthY, lengthZ
	   Output:	- points of the mesh are translated
	   	   			
	6. Curvature 2D
	   Calculate the curvature for each point by summing up all the angles of 
	   the triangles around the given point and substract this sum of 2*PI in 
	   dimension 2
	   The function creates the file sol in the current working directory 
	   Input:	- mesh mesh
	   Output:  - Curvature for each point in mesh->sol
	   
	7. Curvature 3D
	   Calculate the curvature for each point by summing up all the angles of 
	   the triangles around the given point and substract this sum of 2*PI in 
	   dimension 2
	   The function creates the file sol in the current working directory 
	   Input:	- mesh mesh
	   Output:  - Curvature for each point in mesh->sol
	   
	8. Bucket
	   Subdivise the mesh in several sub-domaines of the same size.
	   This allows us to decrease highly the calculation time and complexity
	   of certain calculations 
	   Steps:  - Creation of the bucket structure
	   		   - Fill up bucket
	   		   - Using the bucket
	   			
	9. Distance point to triangle
	   Calculate the distance between a point and a triangle by using
	   the barycentric coordinates
	   Input:	- mesh mesh
	   		    - triangle
	   		    - point
	   Output:  - distance
	   
	10. Hash function
	    Creation of the adjacent relations of a mesh 	   
	    ( Several fonctions do this job, but the user just uses the
	     function setAdj)
	    Input:	- mesh mesh
	    Output:	- adjacent relations
	    
	11. Normales
		Calculate the normales of the triangles of a mesh with the help 
		of a cross product
		Input:	- mesh mesh
		Output:	- normales are saved in mesh->nn and in the file 
		          "normales.mesh"
				  
	12. Distance point to mesh via bucket
		Calculate the distance between a point and a mesh with the help of 
		the bucket
		Input:	- mesh mesh
				- point
				- VertToTria (Array which has as index the vertices and as 
				  value the corresponding triangle)
				- bucket
		Output:	- Distance
		
	13. Ball (check that the given point does exist in the mesh!)
		Calculate the ball of triangles around a given point with bruteforce 
		and the method of adjacences
		Input:	- mesh mesh
				- start triangle
				- start point
				- pointer to list in which the found triangles are saved 
		Output:	- The found triangles are saved in the given list 
				
	14. distance Hausdorff
		calculate the Hausdorff distance between two given mesh 'meshA', 
		'meshB'
		Input:	- meshs meshA, meshB
		Output:	- Hausdorff distance
	
	
	

************************ Visualisation ***********************************


To visualise the files .mesh:

1. Download and install "medit 3.0" in the rubric "SCIENTIFIC VISUALIZATION" 
   of the web page: "http://www.ann.jussieu.fr/frey/software.html"

2. Start the program with the command "medit-linux" and the name of file
   which should be shown (with path if necessary)
   
   Exemple: "medit-linux ../mesh/cube.mesh"

