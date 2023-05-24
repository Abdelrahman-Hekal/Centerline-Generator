import sys
import os
sys.path.append(os.path.dirname(__file__))
import IO
import math
import copy
import sys
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d
#import numpy as np

point_nfaces = {}
point_face_ids = {}
faces_edges = {}

def process_fluent_mesh(path, name='inner'):

    """
    inputs:
    1. pathectory for fluent ascii mesh file
    2. boundary name to be extracted from the mesh (optional), the default case is the name that contains the word 'inner'

    outputs:
    1. VTK output for the extracted boundary
    2. points and faces of the extracted boundary
    """
    boundaries = []
    #reading fluent input file
    points, cells = IO.read_ASCII_fluent_mesh(path)
    points = points[0]

    boundary = []
    boundary_points = []
    points_updated_id = {}
    inner_boundary_id = 0
    points_id = {}
    boundary_found = False
    # searching the inner boundary id
    for key, value in IO.boundaries.items():
        if name.lower() == key.lower():
            inner_boundary_id = value
            boundary_found = True
            break
    if not boundary_found:
        print('The given boundary doesn\'t exist in the mesh')
        exit(1)
    # extracting the desired boundary
    for i, id in enumerate(IO.faces_ids):
        if id == inner_boundary_id:
            for cell in cells[i][1]:   # looping at each face of the inner boundary
                for index, id in enumerate(cell):   
                    if points_id.get(id, 0) == 0:   # first time to process point
                        boundary_points.append(points[id])  # adding the point with a new id starts from zero
                        cell[index] = len(boundary_points) - 1  # updating the point id for the face
                        points_id[id] = 1                       # indicates that the point is already processed
                        points_updated_id[id] = len(boundary_points) - 1    # mapping from the old to the new ids
                    else:
                        cell[index] = points_updated_id[id] # updating the point id for the face
                        points_id[id] += 1  # counter for detecting boundary points

            # considering the faces of the inner boundary only
            boundary.append(cells[i])
            
    # updating the cells and points variables with the inner boundary data only
    cells.clear()
    points.clear()
    cells = boundary[:]
    #points = boundary_points[:]

    # output vtk file for testing
    filename = path + '_' + name + '_geo''.vtk'
    IO.unstructured_grid_VTK_writer(filename, boundary_points, cells)

    faces = {}
    i = 0
    for j in range(len(cells)):
        for face in cells[j][1]:
            faces[i] = face
            i += 1

    nodes = {}
    for i, node in enumerate(boundary_points):
            nodes[i] = node

    return nodes, faces
   
def centerline_generation(points, faces, id1=0, id2=1):

    """
    inputs:
    1. points and faces of a boundary

    outputs:
    1. points and radii of the centerline for the given boundary
    """
    ports_edges = []
    ports_centers = []
    boundary_edges = []
    ports = [[],[]]
    prev_ports = [[],[]]
    centerline_tmp = []
    radii_tmp = []
    freeze = []
    boundary_faces = []
    prev_boundary_faces = []
    done = False
    
    boundary_edges = get_boundaries(points, faces, freeze)
    nports, ports_edges = ports_detection(boundary_edges)

    if nports > 2:
        ports_centers = get_port_center1(points, ports_edges)
        ports, freeze, ports_centers = get_ports_from_user(ports_centers, ports_edges, centerline_tmp, id1, id2)
        get_radii(ports, points, radii_tmp, ports_centers)
    elif nports < 2:
        print('The given boundary must have at least two ports, centerline generation failed')
        exit(1)
    else:
        ports_centers = get_port_center1(points, ports_edges)
        for center in ports_centers:
            centerline_tmp.append(center)

        for i, port in enumerate(ports_edges):
            for edge in port:
                ports[i].append(edge[0])

        get_radii(ports, points, radii_tmp, ports_centers)

    
    while len(faces) > len(prev_ports[0])*2:
        prev_ports, ports, boundary_faces, prev_boundary_faces, done = update_ports(prev_ports, ports, freeze, boundary_faces, done, nports)

        if len(ports[0]) == 0 or len(ports[1]) == 0 or done:
            break
        ports_centers = get_port_center2(points, ports)
        for center in ports_centers:
            centerline_tmp.append(center)

        get_radii(ports, points, radii_tmp, ports_centers)
        update_mesh(faces, points, ports, prev_boundary_faces)


    centerline = [0]*len(centerline_tmp)
    radii = [0]*len(centerline_tmp)

    # ordering centerline points:
    index = 0
    for i in range(len(centerline_tmp)):
        if i % 2 == 0:
            centerline[len(centerline_tmp) - index - 1] = centerline_tmp[i]     
            radii[len(radii_tmp) - index - 1] = radii_tmp[i] 
        else:
            centerline[index] = centerline_tmp[i]
            radii[index] = radii_tmp[i]
            index += 1

    final_line = []
    final_radii = []
    #smoothing centerline
    #######################################################################
    x_coo = []
    y_coo = []
    z_coo = []

    for point in centerline:
        x_coo.append(point[0])
        y_coo.append(point[1])
        z_coo.append(point[2])
     
    sigma = 3
    x_fine, y_fine, z_fine, r_fine = gaussian_filter1d([x_coo, y_coo, z_coo, radii], sigma=sigma, mode='nearest')
  
    if len(centerline) > 50:
        for i in range(2):
            x_fine[0] = x_coo[0]
            x_fine[-1] = x_coo[-1]
            y_fine[0] = y_coo[0]
            y_fine[-1] = y_coo[-1]
            z_fine[0] = z_coo[0]
            z_fine[-1] = z_coo[-1]
            x_fine, y_fine, z_fine, r_fine = gaussian_filter1d([x_fine, y_fine, z_fine, r_fine], sigma=sigma, mode='nearest')

    n = len(x_fine)
    point = []
    for i in range(n):
        point = [x_fine[i], y_fine[i], z_fine[i]]
        final_line.append(point)
        final_radii.append(r_fine[i])
   
    final_line[0] = centerline[0]
    final_line[-1] = centerline[-1]
    ########################################################################
    #smoothed_centerline = Gaussian_smooth(centerline,5)
    #final_line.append(centerline[0])
    #final_line += smoothed_centerline
    #final_line.append(centerline[-1])

    #updating the centerline radii to match the smoothed line
    #j = 0
    #final_radii = []
    #x = len(final_line)
    #for i in range(1, x + 1):
    #    j = (i - 1)/(x - 1)
    #    final_radii.append(((1-j)*radii_tmp[-1] + j*radii_tmp[0]))
              
    return final_line, final_radii

def get_radii(ports, points, radii_tmp, ports_centers):

    for i, port in enumerate(ports):
        r = 1000000  
        for id in port:
            distance = (((ports_centers[i][0] - points[id][0])**2 + (ports_centers[i][1] - points[id][1])**2 + (ports_centers[i][2] - points[id][2])**2)**0.5)
            r = min(distance,r)
        radii_tmp.append(r)


def get_distance(ports_centers, points):

    distance = (((ports_centers[0][0] - ports_centers[1][0])**2 + (ports_centers[0][1] - ports_centers[1][1])**2 + (ports_centers[0][2] - ports_centers[1][2])**2)**0.5)

    return distance

def get_port_center1(points, ports_edges):

    ports_centers = []
    for port in ports_edges:
        center = [0, 0, 0]
        for edge in port:
            center[0] += points[edge[0]][0]
            center[1] += points[edge[0]][1]
            center[2] += points[edge[0]][2]

        center[0] /= float(len(port))
        center[1] /= float(len(port))
        center[2] /= float(len(port))

        ports_centers.append(center)

    return ports_centers

def get_port_center2(points, ports):

    ports_centers = []
    for port in ports:
        center = [0, 0, 0]
        for id in port:
            center[0] += points[id][0]
            center[1] += points[id][1]
            center[2] += points[id][2]

        center[0] /= float(len(port))
        center[1] /= float(len(port))
        center[2] /= float(len(port))

        ports_centers.append(center)

    return ports_centers

def get_boundaries(points, faces, freeze):

    boundary_faces = []
    boundary_edges = []
    edges = []
    edge_face_ids = {}

    # Assessment of faces edges
    for i, face in faces.items():
        npoints = len(face)
        edges.clear()
        for j in range(npoints):
            if j == (npoints - 1):
                edge = (face[j], face[0])
            else:
                edge = (face[j], face[j+1])

            if faces_edges.get(i, 0) == 0:
                faces_edges[i] = [edge]
            else:
                faces_edges[i].append(edge)

            point_nfaces[face[j]] = point_nfaces.get(face[j], 0) +  1
            if point_face_ids.get(face[j], 0) == 0:
                point_face_ids[face[j]] = [i]
            else:
                point_face_ids[face[j]].append(i)

            edge2 = (edge[1], edge[0])
            if edge_face_ids.get(edge2, 0) == 0:
                if edge_face_ids.get(edge, 0) == 0:
                    edge_face_ids[edge] = [i]
                else:
                    edge_face_ids[edge].append(i)
            else:
                edge_face_ids[edge2].append(i)

    # Assessment of boundary edges and faces
    for edge, face_ids in edge_face_ids.items():
        if len(face_ids) == 1:
            if face_ids[0] not in boundary_faces:
                boundary_faces.append(face_ids[0])
            boundary_edges.append(edge)
            for id in edge:
                for face in point_face_ids[id]:
                    if face not in boundary_faces:
                        boundary_faces.append(face)

    return boundary_edges

def ports_detection(boundary_edges):

    if len(boundary_edges) < 6:
        print('Each Port Must Have At Least 3 Edges, Centerline Generation Failed')
        exit(1)

    nports = 0
    port_edges = []
    ports_edges = []
    boundary_edges_tmp = boundary_edges.copy()
    while len(boundary_edges_tmp) > 2:
        nports += 1
        port_found = False
        port_edges.clear()
        port_edges.append(boundary_edges_tmp[0])
        edge1 = port_edges[0]
        prev_edge = ()
        # port points detection
        while not port_found:
            for edge in boundary_edges_tmp:
                if edge == edge1: continue
                if edge[0] in edge1 or edge[1] in edge1: 
                # condition for closing the port contour
                    if edge == port_edges[0] and len(port_edges) > 2:
                        port_found = True
                        for edge in port_edges:
                            boundary_edges_tmp.remove(edge)
                        ports_edges.append(port_edges[:])
                        break
                    elif edge not in port_edges:
                        done = False
                        port_edges.append(edge)
                        edge1 = edge
                        continue
                
    return nports, ports_edges


def get_nearst_id(point1, list, points):

    if len(list) == 0:
        print('Centerline Generation Failed!')
        exit(1)
    clearance = {}
    for id in list:
        point2 = points[id]
        dis = ((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2 + (point2[2] - point1[2])**2)**0.5
        if dis > 0:
            clearance[id] = dis
                        
    min_dis = min(clearance.values())

    for id in clearance.keys():
        if clearance[id] == min_dis:
            return id

def Gaussian_smooth(list1, degree=5):
    """
    inputs:
    1. list of points coordinates to be smoothed
    2. degree of smoothness

    output:
    1. list of smoothed points coordinates
    """
    window=(degree*2-1)
    weight=[1.0]*window
    weightGauss=[]  
    for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(math.exp((4*(frac))**2))  
         weightGauss.append(gauss)  

    weight = [a * b for a, b in zip(weightGauss, weight)]
    total = sum(weight)
    smoothed=[0.0]*(len(list1)-window)

    for i in range(len(smoothed)):
        smoothing=[0.0, 0.0, 0.0]
        for e,w in zip(list1[i:i+window],weight):
            for j in range(len(smoothing)):
                smoothing[j] += w*e[j]
        for k in range(len(smoothing)):
            smoothing[k] /= total
        smoothed[i]=smoothing


    return list(smoothed)
 
def update_mesh(faces, points, ports, prev_boundary_faces):

    for port in ports:
        for id in port:
            face_ids = point_face_ids[id].copy()
            for face in face_ids:
                if face in prev_boundary_faces:
                    values = point_face_ids[id]
                    values.remove(face)
                    if face in faces:
                        faces.pop(face)


def update_ports(prev_ports, ports, freeze, boundary_faces, done, nports):

    prev_ports = copy.deepcopy(ports)
    prev_boundary_faces = copy.deepcopy(boundary_faces)
    boundary_faces.clear()
    ports[0].clear()
    ports[1].clear()
    for i, port in enumerate(prev_ports):
        null_point = 0
        for id in port:
            if len(point_face_ids[id]) == 0:
                null_point += 1
                if (nports == 2 and null_point > 8) or (nports > 2 and null_point > 1):
                    done = True
            for face in point_face_ids[id]:
                for edge in faces_edges[face]:
                    for elem in freeze:
                        if edge in elem: continue
                    if edge[0] not in prev_ports[i] and edge[0] not in ports[i]:
                        ports[i].append(edge[0])
                        for face2 in point_face_ids[edge[0]]:
                            if face2 not in prev_boundary_faces and face2 not in boundary_faces:
                                boundary_faces.append(face2)

    return prev_ports, ports, boundary_faces, prev_boundary_faces, done

def get_ports_from_user(ports_centers, ports_edges, centerline_tmp, id1, id2):
    
    try:
        id1 = int(id1)
        id2 = int(id2)
    except:
        print('The parsed points IDs are invalid!')
        exit(2)

    if id1 > len(ports_centers) - 1 or id2 > len(ports_centers) - 1:
        print('The parsed points IDs are invalid!')
        exit(2)
    selected_ports = []
    # inlet and outlet ports specified by the user
    selected_ports.append(ports_centers[id1])
    selected_ports.append(ports_centers[id2])

    ports_centers_tmp = []
    ports = [[],[]]
    freeze = []
    j = 0
    for i, center in enumerate(ports_centers):
        matched = False
        if i == id1 or i == id2:
            matched = True
            for edge in ports_edges[i]:
                ports[j].append(edge[0])
            j += 1  
            centerline_tmp.append(center)
            continue
        if not matched:
            if ports_edges[i] not in freeze:
                freeze.append(ports_edges[i])     

    ports_centers.clear()
    ports_centers = copy.deepcopy(centerline_tmp)

    if len(ports_centers) != 2:
        print('Parsed ports centers are invlid')
        exit(1)

    return ports, freeze, ports_centers

def detect_ports(path, boundary):

    ports_edges = []
    ports_centers = []
    boundary_edges = []
    freeze = []

    points, faces = process_fluent_mesh(path, boundary)
    boundary_edges = get_boundaries(points, faces, freeze)
    nports, ports_edges = ports_detection( boundary_edges)
    ports_centers = get_port_center1(points, ports_edges)

    if nports != len(ports_centers):
        print('Ports detection failed')
        exit(2)

    print(len(ports_centers))
    print(ports_centers)
    return ports_centers
            
def generate_centerline(path, boundary='inner', id1=0, id2=1):
    """
    inputs:
    1. pathectory of ASCII fluent mesh file
    2. boundary label to be used in the centerline generation
    3. port1: an optional argument for the inlet port center in case of multi-port geometries
    3. port2: an optional argument for the outlet port center in case of multi-port geometries

    outputs:
    1. VTK output for the extracted boundary
    2. VTK output for the centerline of the given boundary
    """
    filename = path.split("\\")[-1]
    points, faces = process_fluent_mesh(path, boundary)
    centerline, radii = centerline_generation(points, faces, id1, id2)
    IO.polyline_VTK_writer(centerline, radii, path, filename, boundary)

##############################################################################################################
# the main functions to be called:
path = 'D:/Hekal/VS projects/centerline generation/Centerline Generation/test cases/1658063426280-6.PPAFTDOC.stl.msh'
boundary = 'Inner_C_2'
id1 = '0'
id2 = 1

#detect_ports(path, boundary)
generate_centerline(path, boundary, id1, id2)

#if __name__ == '__main__':
#   # print (len(sys.argv))
#   params = []
#   if len(sys.argv) > 2:
#      for x in range(2, len(sys.argv)):
#         params.append(sys.argv[x])
#   if len(params) > 0:
#      globals()[sys.argv[1]](*params)
#   else:
#      globals()[sys.argv[1]]()
#   # input()




