******************* Démarrer le programme ****************************

Etapes:

A. Allez dans le répertoire ..Polytech_maillage/prog

B. Compilez le programme avec la commande "make"

C. Lancez le programme avec la commande "./main"

	Maintenant vous pouvez choisir entre les cas de test suivants:
	
	1. Rotation 2D
	   rotation via l'angle 'angle' dans dimension 2
	   Input:	- maillage mesh
	   			- angle
	   Output:	- Points enregistrés dans le mesh sont tournés
	   
	2. Rotation 3D
	   rotation via les angles angleX, angleY, angleZ dans dimension 3
	   Input:	- maillage mesh
	   			- angle
	   Output:	- Points enregistré dans le mesh sont tournés
	   			
	3. Superposition
	   Fusion tous les informations (points, edge, triangles,...) de deux
	   maillage Mesh1 et Mesh2 dans un nouveau maillage Mesh_final
	   Input:	- maillage Mesh1, Mesh2, Mesh_final)
	   Output:  - tous les informations sont enregistrées dans Mesh_final
	   
	4. Translation 2D
	   déplace chaque point du maillage par les valeur lengthX, lengthY
	   en dimension 2
	   Input: 	- maillage mesh
	   			- lengthX, lengthY
	   Output:	- Points enregistrés dans le mesh sont déplacés
	   			
	5. Translation 3D
	   déplace chaque point dans le maillage par les valeur lengthX, lengthY,
	   lengthZ en dimension 3
	   Input: 	- maillage mesh
	   			- lengthX, lengthY, lengthZ
	   Output:	- Points enregistrés dans le mesh sont déplacés 	   			
	6. Curvature 2D
	   Calcule la courbure pour chaque point en sommant tous les angles des 
	   triangles autour d'un point et soustrait cette somme de 2*PI en 
	   dimension 2
	   La fonction crée le fichier sol dans répértoire courant 
	   Input:	- maillage mesh
	   Output:  - Courbure par point dans mesh->sol
	   
	7. Curvature 3D
	   Calcule la courbure pour chaque point en sommant tous les angles des 
	   triangles autour d'un point et soustrait cette somme de 2*PI en 
	   dimension 3
	   La fonction crée le fichier sol dans répértoire courant 
	   Input:	- maillage mesh
	   Output:  - Courbure par point dans mesh->sol	
	   
	8. Bucket
	   Subdivise le maillage en plusieurs sous-domaines de même taille.
	   Ca nous permet de diminuer grandement le temps et la lourdeur de 
	   certains calculs
	   Etapes: - creation de la structure bucket
	   		   - remplissage du bucket
	   		   - utilisation du bucket
	   			
	9. Distance point to triangle
	   calcule la distance entre un point et un triangle en utilisant les 
	   coordonnées barycentriques
	   Input:	- maillage mesh
	   		    - triangle
	   		    - point
	   Output:  - distance
	   
	10. Hash function
	    Creation des relations d'adjacence à partir d'un maillage 	   
	    ( Plusieurs fonctions réalisées, seule la fonction setAdj est 
	     utilisée)
	    Input:	- maillage mesh
	    Output:	- relations adjacence
	    
	11. Normales
		Calcule les normales des triangles d'un maillage avec l'aide d'un 
		produit vectoriel
		Input:	- maillage mesh
		Output:	- normales sont enregistrées dans mesh->nn et dans un 
				  fichier "normales.mesh"
				  
	12. Distance point to mesh via bucket
		Calcule la distance entre un point au maillage via le bucket
		Input:	- maillage mesh
				- point
				- VertToTria (Tableau qui a comme indice les vertices et qui 
				  a associé un des triangles correspondant)
				- bucket
		Output:	- Distance
		
	13. Ball (check que le point donné exist dans le maillage!)
		Calcule boule des triangles autour un point donné avec bruteforce et 
		avec la méthode d'adjacence
		Input:	- maillage mesh
				- start triangle
				- start point
				- pointer to liste dans laquelle les triangles trouvés sont 
				  enregistré
		Output:	- Les triangles trouvés sont enregistrés dans la liste donnée 
				
	14. distance Hausdorff
		calcule la distance de Hausdorff entre deux maillage donnés meshA, 
		meshB
		Input:	- maillages meshA, meshB
		Output:	- distance de Hausdorff
	
	
	

************************ Visualiser ***********************************


Pour visualiser les fichier .mesh:

1. Téléchargez "medit 3.0" dans la rubrique "SCIENTIFIC VISUALIZATION" de 
   la site web: http://www.ann.jussieu.fr/frey/software.html

2. Lancez le programme avec la commande "medit-linux" et le nom de votre 
   fichier (avec chemin d'accès si nécessaire)
   
   Exemple: "medit-linux ../mesh/cube.mesh"

