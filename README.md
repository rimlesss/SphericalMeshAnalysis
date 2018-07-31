# SphericalMeshAnalysis
Spherical parameterization of triangle meshes, with a goal of minimizing area distortion

first goal: initial parametrization
- read mesh
- check for genus
- (optional) preprocessing of mesh
- fix poles from mesh
Now following Gotsman 2003 and Saba 2005 
- initialize a symmetric Laplacian matrix 

(Optional) Alternative method based on progressive mesh, following Praun 2003