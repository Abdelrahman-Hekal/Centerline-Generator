import re

# dict for storing the mesh wall boundaries and ids
boundaries = {}
# list for storing the associated boundary id for the processd boundary
faces_ids = []

"""
I/O for Ansys's mesh, VTK centerlines and surfaces 

"""

def _skip_to(f, char):
    c = None
    while c != char:
        c = f.read(1).decode("utf-8")
    return

def _skip_close(f, num_open_brackets):
    while num_open_brackets > 0:
        char = f.read(1).decode("utf-8")
        if char == "(":
            num_open_brackets += 1
        elif char == ")":
            num_open_brackets -= 1
    return

def _read_points(f, line, first_point_index_overall, last_point_index):
    # If the line is self-contained, it is merely a declaration
    # of the total number of points.
    if line.count("(") == line.count(")"):
        return None, None, None

    # (3010 (zone-id first-index last-index type ND)
    out = re.match("\\s*\\(\\s*(|20|30)10\\s*\\(([^\\)]*)\\).*", line)
    a = [int(num, 16) for num in out.group(2).split()]
    a.append(3) #update
    if len(a) <= 4:
        raise ReadError()

    first_point_index = a[1]
    # store the very first point index
    if first_point_index_overall is None:
        first_point_index_overall = first_point_index
    # make sure that point arrays are subsequent
    if last_point_index is not None:
        if last_point_index + 1 != first_point_index:
            raise ReadError()
    last_point_index = a[2]
    num_points = last_point_index - first_point_index + 1
    dim = a[4]

    # Skip ahead to the byte that opens the data block (might
    # be the current line already).
    last_char = line.strip()[-1]
    while last_char != "(":
        last_char = f.read(1).decode("utf-8")

    if out.group(1) == "":
        # ASCII data
        pts = [0]*num_points
        for i in range(len(pts)):
            pts[i] = [0]
            for j in range(dim):
                if j == 0: continue
                pts[i].append(0)

        for k in range(num_points):
            # skip ahead to the first line with data
            line = ""
            while line.strip() == "":
                line = f.readline().decode("utf-8")
            dat = line.split()
            if len(dat) != dim:
                raise ReadError()
            for d in range(dim):
                pts[k][d] = float(dat[d])

    # make sure that the data set is properly closed
    _skip_close(f, 2)
    return pts, first_point_index_overall, last_point_index

def _read_cells(f, line):
    # If the line is self-contained, it is merely a declaration of the total number of
    # points.
    if line.count("(") == line.count(")"):
        return None, None

    out = re.match("\\s*\\(\\s*(|20|30)12\\s*\\(([^\\)]+)\\).*", line)
    a = [int(num, 16) for num in out.group(2).split()]
    if len(a) <= 4:
        raise ReadError()
    first_index = a[1]
    last_index = a[2]
    num_cells = last_index - first_index + 1
    zone_type = a[3]
    element_type = a[4]

    if zone_type == 0:
        # dead zone
        return None, None

    key, num_nodes_per_cell = {
        0: ("mixed", None),
        1: ("triangle", 3),
        2: ("tetra", 4),
        3: ("quad", 4),
        4: ("hexahedron", 8),
        5: ("pyra", 5),
        6: ("wedge", 6),
    }[element_type]

    # Skip to the opening `(` and make sure that there's no non-whitespace character
    # between the last closing bracket and the `(`.
    if line.strip()[-1] != "(":
        c = None
        while True:
            c = f.read(1).decode("utf-8")
            if c == "(":
                break
            if not re.match("\\s", c):
                # Found a non-whitespace character before `(`.
                # Assume this is just a declaration line then and
                # skip to the closing bracket.
                _skip_to(f, ")")
                return None, None

    if key == "mixed":
        data = None
    else:
        # read cell data
        if out.group(1) == "":
            # ASCII cells
            data = [0]*num_nodes_per_cell
            for i in range(len(data)):
                data[i] = [0]
                for j in range(num_cells):
                    if j == 1: continue
                    data[i].append(0)

            for k in range(num_cells):
                line = f.readline().decode("utf-8")
                dat = line.split()
                if len(dat) != num_nodes_per_cell:
                    raise ReadError()
                data[k] = [int(d, 16) for d in dat]

    # make sure that the data set is properly closed
    _skip_close(f, 2)
    return key, data

def _read_faces(f, line):
    # faces
    # (13 (zone-id first-index last-index type element-type))
    # If the line is self-contained, it is merely a declaration of
    # the total number of points.
    id = 0
    if line.count("(") == line.count(")"):
        return {}

    out = re.match("\\s*\\(\\s*(|20|30)13\\s*\\(([^\\)]+)\\).*", line)
    a = [int(num, 16) for num in out.group(2).split()]
    id = int(line[4:][line[4:].index('(') + 1:line[4:].index(' ')], 16)
    #try:
    #    id = int(out.group(0)[5:8], 16)
    #except ValueError:
    #    id = int(out.group(0)[5:7], 16)

    
    if len(a) <= 4:
        raise ReadError()
    first_index = a[1]
    last_index = a[2]
    num_cells = last_index - first_index + 1
    element_type = a[4]

    element_type_to_key_num_nodes = {
        0: ("mixed", None),
        2: ("line", 2),
        3: ("triangle", 3),
        4: ("quad", 4),
    }

    key, num_nodes_per_cell = element_type_to_key_num_nodes[element_type]
    # Skip ahead to the line that opens the data block (might be
    # the current line already).
    if line.strip()[-1] != "(":
        _skip_to(f, "(")

    data = {}
    if out.group(1) == "":
        # ASCII
        if key == "mixed":
            # > If the face zone is of mixed type (element-type = > 0), the body of the
            # > section will include the face type and will appear as follows
            # > type v0 v1 v2 c0 c1
            for k in range(num_cells):
                line = ""
                while line.strip() == "":
                    line = f.readline().decode("utf-8")
                dat = line.split()
                type_index = int(dat[0], 16)
                if type_index == 0:
                    raise ReadError()
                type_string, num_nodes_per_cell = element_type_to_key_num_nodes[
                    type_index
                ]
                if len(dat) != num_nodes_per_cell + 3:
                    raise ReadError()

                if type_string not in data:
                    data[type_string] = []

                data[type_string].append(
                    [int(d, 16) for d in dat[1 : num_nodes_per_cell + 1]]
                )

            data = {key: data[key] for key in data}

        else:
            # read cell data
            data = [0]*num_cells
            for i in range(len(data)):
                data[i] = [0]
                for j in range(num_nodes_per_cell):
                    if j == 1: continue
                    data[i].append(0)

            for k in range(num_cells):
                line = f.readline().decode("utf-8")
                dat = line.split()
                # The body of a regular face section contains the grid connectivity, and
                # each line appears as follows:
                #   n0 n1 n2 cr cl
                # where n* are the defining nodes (vertices) of the face, and c* are the
                # adjacent cells.
                if len(dat) != num_nodes_per_cell + 2:
                    raise ReadError()
                data[k] = [int(d, 16) for d in dat[:num_nodes_per_cell]]
            data = {key: data}

    for i in range(len(data.keys())):
        faces_ids.append(id)
    # make sure that the data set is properly closed
    _skip_close(f, 2)

    return data

def read_ASCII_fluent_mesh(filename):  
    # Initialize the data optional data fields
    field_data = {}
    cell_data = {}
    point_data = {}

    points = []
    cells = []

    first_point_index_overall = None
    last_point_index = None

    # read file in binary mode since some data might be binary
    with open (filename, "rb") as f:
        while True:
            line = f.readline().decode("utf-8")
            if not line:
                break

            if line.strip() == "":
                continue

            # expect the line to have the form
            #  (<index> [...]
            out = re.match("\\s*\\(\\s*([0-9]+).*", line)
            if not out:
                print('unsupported mesh format')
                exit(1)
            index = out.group(1)

            if index == "0":
                # Comment.
                _skip_close(f, line.count("(") - line.count(")"))
            elif index == "1":
                # header
                # (1 "<text>")
                _skip_close(f, line.count("(") - line.count(")"))
            elif index == "2":
                # dimensionality
                # (2 3)
                _skip_close(f, line.count("(") - line.count(")"))
            elif re.match("(|20|30)10", index):
                # points
                pts, first_point_index_overall, last_point_index = _read_points(
                    f, line, first_point_index_overall, last_point_index
                )

                if pts is not None:
                    points.append(pts)

            elif re.match("(|20|30)12", index):
                # cells
                # (2012 (zone-id first-index last-index type element-type))
                key, data = _read_cells(f, line)
                if data is not None:
                    cells.append((key, data))

            elif re.match("(|20|30)13", index):
                data = _read_faces(f, line)

                for key in data:
                    cells.append((key, data[key]))

            elif index == "39":
                obj = re.match("\\(39 \\([0-9]+ ([\\S]+) ([\\S]+)\\)\\(\\)\\)", line)
                _skip_close(f, line.count("(") - line.count(")"))
                if obj:
                    boundary_id = re.findall(r'\d+', line)[1]
                    boundaries[obj.group(2)] =  int(boundary_id)
                    #try:
                    #    boundaries[obj.group(2)] =  int(obj.group(0)[5:9])
                    #except ValueError:
                    #    try:
                    #        boundaries[obj.group(2)] =  int(obj.group(0)[5:8])
                    #    except ValueError:
                    #        try:
                    #            boundaries[obj.group(2)] =  int(obj.group(0)[5:7])
                    #        except ValueError:
                    #            boundaries[obj.group(2)] =  int(obj.group(0)[5:6])

            elif index == "45":
                # (45 (2 fluid solid)())
                obj = re.match("\\(45 \\([0-9]+ ([\\S]+) ([\\S]+)\\)\\(\\)\\)", line)
                
                if obj:
                    boundary_id = re.findall(r'\d+', line)[1]
                    boundaries[obj.group(2)] =  int(boundary_id)
                    #try:
                    #    boundaries[obj.group(2)] =  int(obj.group(0)[5:9])
                    #except ValueError:
                    #    try:
                    #        boundaries[obj.group(2)] =  int(obj.group(0)[5:8])
                    #    except ValueError:
                    #        try:
                    #            boundaries[obj.group(2)] =  int(obj.group(0)[5:7])
                    #        except ValueError:
                    #            boundaries[obj.group(2)] =  int(obj.group(0)[5:6])
            else:
                _skip_close(f, line.count("(") - line.count(")"))

    # Gauge the cells with the first point_index.
    for list in cells:
        for cell in (list[1]):
            for j in range(len(cell)): 
                cell[j] -= first_point_index_overall

    return points, cells

def polyline_VTK_writer(centerline, radii, dir, filename, boundary):
    if not centerline:
        print('Empty line error!')
        exit(1)
    label = '_' + boundary + '_centerline'
    file = dir + label + '.vtk'
    #writing file header
    with open (file, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('vtk output\n')
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS {} float\n'.format(len(centerline)))
    
        #writing points
        value = 0
        i = 0
        j = 0
        for point in centerline:
            j += 1
            if  j == 1:
                end_point = '{:.2f} {:.2f} {:.2f} '.format(point[0], point[1], point[2])
            elif j == len(centerline):
                f.write(end_point)
                value += 1
                if value % 3 == 0:
                    f.write('\n')
                string = '{:.2f} {:.2f} {:.2f} '.format(point[0], point[1], point[2])
                f.write(string)
                value += 1
                if value % 3 == 0:
                    f.write('\n')
            else:
                string = '{:.2f} {:.2f} {:.2f} '.format(point[0], point[1], point[2])
                f.write(string)
                value += 1
                if value % 3 == 0:
                    f.write('\n')
        
        f.write('\n')

        #writing points order
        f.write('LINES 1 {}\n'.format(len(centerline) + 1))
        for i in range(len(centerline), -1, -1):
            if i == (len(centerline) - 2):
                string2 = '{} '.format(i)
            else:
                string = '{} '.format(i)
                f.write(string)
        f.write(string2)

        #writing radii
        f.write('\n')
        f.write('POINT_DATA {}\n'.format(len(centerline)))
        f.write('FIELD FieldData 1\n')
        f.write('MaximumInscribedSphereRadius 1 {} double\n'.format(len(centerline)))

        value = 0
        j = 0
        for i in range(len(centerline)):
            r= radii[i]
            j += 1
            if  j == 1:
                end_point = '{:.2f} '.format(r)
            elif j == len(centerline):
                f.write(end_point)
                value += 1
                if value % 9 == 0:
                    f.write('\n')
                string = '{:.2f} '.format(r)
                f.write(string)
                value += 1
                if value % 9 == 0:
                    f.write('\n')
            else:
                string = '{:.2f} '.format(r)
                f.write(string)
                value += 1
                if value % 9 == 0:
                    f.write('\n')
        f.write('\n')

def unstructured_grid_VTK_writer(file, boundary_points, cells):

    #writing file header
    with open (file, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('vtk output\n')
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS {} double\n'.format(len(boundary_points)))
    
        #writing points
        value = 0
        for point in boundary_points:
                string = '{:.2f} {:.2f} {:.2f} '.format(point[0], point[1], point[2])
                f.write(string)
                value += 1
                if value % 3 == 0:
                    f.write('\n')
        
        f.write('\n')
        total_ids = 0
        nfaces = 0
        for list in cells:
                total_ids += (len(list[1][0]) + 1)*len(list[1])
                nfaces += len(list[1])

        
        #writing faces
        f.write('CELLS {} {}\n'.format(nfaces, total_ids))

        for list in cells:
            for face in list[1]:
                string = ''
                string += str(len(face)) + ' '
                for id in face:
                    string += str(id) + ' '
                f.write(string)
                f.write('\n')
        f.write('\n')

        f.write('CELL_TYPES {}\n'.format(nfaces))
        for list in cells:
            for face in list[1]:
                if len(face) == 3:
                    f.write('{}\n'.format(5))
                elif len(face) == 4:
                    f.write('{}\n'.format(9))
                else:
                    f.write('{}\n'.format(7))
