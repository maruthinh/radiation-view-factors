
import numpy as np 
from numba.experimental import jitclass # import the decorator
import numba.types as nbt  # import the types
import pyvista as pv 
import time

spec = [("surf1_vertices", nbt.f8[:,:]), ("surf2_vertices", nbt.f8[:,:]), \
        ("surf1_normal", nbt.f8[:]), ("surf2_normal", nbt.f8[:])]

@jitclass(spec)
class ViewFactor(object):
    """To compute view factors between two rectangular surfaces"""
    def __init__(self, surf1_vertices, surf1_normal, surf2_vertices, surf2_normal):
        self.surf1_vertices = surf1_vertices
        self.surf2_vertices = surf2_vertices
        self.surf1_normal = surf1_normal
        self.surf2_normal = surf2_normal
        
    def random_point_on_rectangle(self):
        """Compute random point on a rectangle surface"""
        v1, v2, v3, v4 = self.surf1_vertices
        xi1 = np.random.rand()
        xi2 = np.random.rand()
        p1 = v1 + xi1*(v2-v1)
        return  p1 + xi2*(v4-v1)
    
    def random_direction(self):
        """Compute random direction on a rectangle surface with hemisphere
        alignment"""
        xi3 = np.random.rand()
        xi4 = np.random.rand()
        theta = np.arcsin(np.sqrt(xi3))
        phi = 2 * np.pi * xi4
        
        direction = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        
        if abs(self.surf1_normal[0]) > abs(self.surf1_normal[1]):
            up = np.array([0.0, 1.0, 0.0])
        else:
            up = np.array([1.0, 0.0, 0.0])
            
        if np.allclose(up, self.surf1_normal):
            return direction
        elif np.allclose(up, -self.surf1_normal):
            return -direction
        else:
            tangent = np.cross(up, self.surf1_normal)
            tangent = tangent / np.linalg.norm(tangent)
            bitangent = np.cross(self.surf1_normal, tangent)
            
            return tangent*direction[0] + bitangent*direction[1] + self.surf1_normal*direction[2]
            
    def check_ray_rectangle_intersection1(self):
        """Check if ray intersects rectangular surface"""
        v1, v2, v3, v4 = self.surf2_vertices 
        rand_dir = self.random_direction()
        rand_pt = self.random_point_on_rectangle()
        denominator = np.dot(rand_dir, self.surf2_normal)
        # Check if ray is parallel to the plane
        if np.abs(denominator) < 1e-10:
            return False
        
        #approaching front 
        if not np.dot(self.surf2_normal, rand_dir) < 0:
            return False 
        
        # compute distance to the plane
        t = np.dot(v1 - rand_pt , self.surf2_normal) / denominator
        if t < 0:
            return False  # Intersection behind the ray origin
        
        # Compute intersection point
        intersection = rand_pt + t * rand_dir
       
        # Check if the intersection point is inside the rectangle using barycentric coordinates or edge tests
        vp = intersection - v1 
        v12 = v2 - v1 
        v14 = v4 - v1
        a = np.dot(vp, v12) / np.dot(v12, v12)
        b = np.dot(vp, v14) / np.dot(v14, v14)
    
        return 0 <= a <= 1 and 0 <= b <= 1    

    def compute_view_factor(self, num_rays):
        """Compute view factors"""
        hits = 0
        for _ in range(num_rays):
            if self.check_ray_rectangle_intersection1():
                hits += 1
        view_factor = hits / num_rays
        return view_factor
            
if __name__ == "__main__": 
    #array of vertices
    vertices = np.array([
        [0.0, 0.0, 0.0], [0.596, 0.0, 0.0], [0.596, 0.0, 0.996], [0.0, 0.0, 0.996],
        [0.0, 1.997, 0.0], [0.596, 1.997, 0.0], [0.596, 1.997, 0.996], [0.0, 1.997, 0.996]
    ])
        
    #connectivity of vertices that define surfaces 
    connectivity = np.array([
        [4, 0, 1, 2, 3],  # bottom
        [4, 4, 5, 6, 7],  # top
        [4, 0, 3, 7, 4],  # left
        [4, 1, 2, 6, 5],  # right
        [4, 3, 2, 6, 7],  # back
        [4, 0, 1, 5, 4]   # front
    ])
    
    #convert to VTK (Pyvista) PolyData object from which we can visualize, and 
    #compute normals, etc 
    connectivity = np.hstack(connectivity)
    surf = pv.PolyData(vertices, connectivity)
    inward_normals=surf.compute_normals(flip_normals=True, inplace=False).cell_normals
    outward_normals=surf.compute_normals(flip_normals=False, inplace=False).cell_normals
    areas=surf.compute_cell_sizes(area=True).get_array('Area')
    
    #visualize surfaces
    surf.plot(
            scalars=areas,
            cpos=[-1, 1, 0.5],
            show_scalar_bar=False,
            show_edges=True,
            line_width=5,
            show_axes=True,
            opacity=0.85,
            screenshot='box.png'
        )
    
    #visualize normals, ray
    p1 = pv.Plotter()
    p1.add_mesh(surf, opacity=0.85, color=True)
    p1.add_axes()
    
    surf1_vertices = np.float64(surf.get_cell(0).points)
    surf1_normal = np.float64(inward_normals[0])
    surf2_vertices = np.float64(surf.get_cell(1).points)
    surf2_normal = np.float64(inward_normals[1])
    vfo = ViewFactor(surf1_vertices, surf1_normal, surf2_vertices, surf2_normal)
    centers = surf.cell_centers().points
    p1.add_arrows(centers[0], inward_normals[0], mag=0.5, color='black', label='surf1_normal')
    p1.add_arrows(centers[1], inward_normals[1], mag=0.5, color='blue', label='surf3_normal')
    p1.add_arrows(vfo.random_point_on_rectangle(), vfo.random_direction(), mag=0.5, color='red', label="rand_dir")
    p1.add_legend()
    
    p1.show()
    
    t1=time.time()
    vf_12 = vfo.compute_view_factor(100000)
    t2=time.time()
    print("view factor from surface 1 to 2 is: \n", vf_12)
    print("time taken is: ", t2-t1)
    
    t1=time.time()
    vf = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            surf1_vertices = np.float64(surf.get_cell(i).points)
            surf1_normal = np.float64(inward_normals[i])
            surf2_vertices = np.float64(surf.get_cell(j).points)
            surf2_normal = np.float64(inward_normals[j])
            vfo = ViewFactor(surf1_vertices, surf1_normal, surf2_vertices, surf2_normal)
            vf[i,j] = vfo.compute_view_factor(100000)
    t2=time.time()
    print("view factor matrix: \n", vf)
    print("time taken is: ", t2-t1)
