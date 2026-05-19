""" Large Displacement Moment Probe for ANSYS Mechanical
    (c) 2026 Yair Preiss

See README and Docs/Manual.md for installation and usage instructions.
"""

import units
from math import sqrt, pi, atan2
from System.Collections.Generic import List
from System import Array


def mag(a):
    """ Return the vector magnitude (norm) of a 3-length list"""  
    
    b = sqrt(a[0]**2.0 + a[1]**2.0 + a[2]**2.0)
    
    return b
    

def create_stm(x, y, z):  
    """ Create a 4x4 spatial transformation matrix 
    
    The spatial transformation matrix maps the effect of a 
    rotation on a vector x, i.e. [M] * x = y. For this implementation,
    translations are not considered.
    
    | x0 x1 x2 0 |
    | y0 y1 y2 0 |
    | z0 z1 z2 0 |
    |  0  0  0 1 |
    
    Args
    ---
    x, y, z: the 3-length rows of the spatial transformation matrix 
        corresponding to the direction of each axis in the global system
    
    Returns
    ---
    4x4 transformation matrix as a linear array in row-major format 
    
    """
    
    matrix = Array.CreateInstance(float,16)
    
    for i in range(3):
        matrix[i] = x[i]
        matrix[i+4] = y[i]
        matrix[i+8] = z[i]
        matrix[i+12] = 0.0
        
    matrix[3] = 0.0
    matrix[7] = 0.0
    matrix[11] = 0
    matrix[15] = 1.0
    
    return matrix        


def cross(a, b):     
    """ Return the cross product of two 3-length lists (vectors)
    """
    
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
         
    return c


def transform(matrix, vector): 
    """ Compute a coordinate transform using the rotation matrix 
    embedded with a transformation matrix
    """
    coordinates = [0.0, 0.0, 0.0]
    
    for i in range(3):
        for j in range(3):
            coordinates[i] += matrix[i][j]*vector[j]
            
    return coordinates


def vsum(a,b):     
    """ Compute the element-wise sum of two vectors, returning a new vector
    """
    
    c = [b[0] + a[0], b[1] + a[1], b[2] + a[2]]
    return c


def avg_disp(points): 
    """ Compute the average displacement within a group of points 
    
    Args 
    ---
    points: list of 3-length lists corresponding to the point coordinates
    
    Returns 
    ---
    list: 3-length vector of x,y,z displacement coordinates 
    """
    
    n = len(points)
    disp = [sum(coord) / n for coord in zip(*points)]
    return disp


def rotation_matrix_from_local_csys(local_csys): 
    """ Return a rotation matrix given a local coordinate system """

    return [local_csys.XAxis, local_csys.YAxis, local_csys.ZAxis]


def get_scale(analysis, is_ld_on):
    """ Get the scalar factors associated with solver and model units """

    reader = analysis.GetResultsData()
    model_units = ExtAPI.DataModel.GeoData.Unit 
    result_locdef = reader.GetResult("LOC_DEF")
    solve_units = result_locdef.GetComponentInfo('X').Unit
    model_scale = units.ConvertToUserUnit(ExtAPI, 1, model_units, "Length")

    if is_ld_on:
        unit_scale = units.ConvertUnit(1, solve_units, model_units,"Length")
    
    else:
        unit_scale = 1.0

    return unit_scale, model_scale


def get_n_elemnodal_forces(mesh, elem_ids, node_ids):
    """ Get the number of elemnodal force results required, by determining which
    nodes are associated with the elements
    """
    count = 0 

    for elem_id in elem_ids:
        element = mesh.ElementById(elem_id)
        for node_id in node_ids:
            if element.NodeIds.IndexOf(node_id) >= 0: 
                # Nodes attached to the element have a positive index
                count += 1
    
    return count


def get_elemnodal_data(
    mesh, 
    result_enfo, 
    result_locdef,
    n_forces, 
    elem_ids, 
    node_ids, 
    is_ld_on, 
    unit_scale
):
    """ 
    """
    node_forces = [[0.0 for _ in range(3)] for _ in range(n_forces)]  
    node_positions = [[0.0 for _ in range(3)] for _ in range(n_forces)]

    count = 0
    for Id in elem_ids:
        element = mesh.ElementById(Id)
        # elementnodal force reactions are reported in a row-major linear array, 
        # with each [Fx1, Fy1, Fz1, Fx2, Fy2, Fz3, ... Fzn], not as individual vectors
        elem_force = result_enfo.GetElementValues(Id) 
        index = List[object]()

        # TODO: this is repeated work, we could just associate node ids 
        # with the element id
        for node_id in node_ids:
            if element.NodeIds.IndexOf(node_id) >= 0:
                index.Add(element.NodeIds.IndexOf(node_id))

        for i in range(len(index)):
            # Assign elementnodal force reaction into node-specific vector
            node_forces[count][0] = elem_force[3*index[i]]
            node_forces[count][1] = elem_force[3*index[i]+1]
            node_forces[count][2] = elem_force[3*index[i]+2]
            if is_ld_on:
                # Assign elementnodal position vectors to pair with respective 
                # force reactions
                node_positions[count] = [
                    x*unit_scale for x in result_locdef.GetNodeValues(
                        element.NodeIds[index[i]]
                        )
                ] 
            else:
                node_positions[count] = [
                    mesh.NodeById(element.NodeIds[index[i]]).X, 
                    mesh.NodeById(element.NodeIds[index[i]]).Y, 
                    mesh.NodeById(element.NodeIds[index[i]]).Z
                ]
            count += 1

    return node_forces, node_positions


def calculate_moment(node_forces, node_positions, centroid, rotation_matrix):
    """ Compute the moment caused by a collection of nodal forces about a 
    specified centroid, in a local coordinate system
    
    Args
    ---
    node_forces: n-length list of [Fx, Fy, Fz] values
        force at each node in global coordinate system
    node_positions: n-length list of [Ux, Uy, Uz] values
        position of each node in global coordinate system (deflected value if 
        NLGEOM=ON)
    centroid: [Cx, Cy, Cz]
        centroid about which to compute the moment 
    rotation_matrix: 3x3 list of lists 
        defines the mapping between the global and local coordinat systems

    Returns 
    ---
    local_moment, r_max : [Mx, My, Mz], float 
        moment in local coordinate system and the maximum distance value, for
        use in plotting
    
    """
    
    n_forces = len(node_forces) 
    assert n_forces == len(node_positions) 
    
    # Compute relative position after deflection (local r vector)
    # Also keep track of largest value for sizing the display vector
    r = [[0.0 for _ in range(3)] for _ in range(n_forces)] 
    global_moment = [0.0, 0.0, 0.0]
    r_max = 0
    for f, force in enumerate(node_forces):
        r[f] = [
            node_positions[f][0] - centroid[0],
            node_positions[f][1] - centroid[1],
            node_positions[f][2] - centroid[2] 
        ]
        if mag(r[f]) > r_max:
            r_max = mag(r[f]) 

        # Sum all M=rxF products in Global orientation (default solver reporting)
        # TODO: logic here could be improved?
        global_moment = vsum(global_moment, cross(r[f], force)) 

    # Transform M into selected CS orientation 
    local_moment = transform(rotation_matrix, global_moment) 
    
    return local_moment, r_max 
    
    

def process_interface(analysis, nodes, local_csys, is_ld_on, unit_scale):
    """ Compute the reaction moment at an interface (i.e. boundary condition)

    Args
    ---
    analysis: ANSYS ACT Analysis object
    nodes: list[int]
        node ids 
    local_csys: ANSYS ACT Coordinate System object
        local coordinate system definition 
    is_ld_on : Bool
        True if NLGEOM is set to ON, False otherwise
    unit_scale: float 
        scale factor for unit system 


    Returns
    ---
    [Mx, My, Mz] in global coordinate system

    """

    reader = analysis.GetResultsData()
    result_locdef = reader.GetResult("LOC_DEF") 
    result_enfo = reader.GetResult("ENFO")
    mesh = analysis.MeshData 
    rotation_matrix = rotation_matrix_from_local_csys(local_csys)
    
    n_nodes = len(nodes)

    # Node positions; this will become a list of 3-length lists
    node_positions = [[0.0 for _ in range(3)] for _ in range(n_nodes)] 

    # Number of element ids is not known a priori 
    elem_ids = []

    # Populate node positions array
    for i, node in enumerate(nodes):
        # Get all elements associated with nodes
        elem_ids += [int(x) for x in mesh.NodeById(node).ConnectedElementIds] 
        if is_ld_on:
            # Get nodal positions at load step
            node_positions[i] = [x * unit_scale for x in result_locdef.GetNodeValues(node)] 
        else:
            # Get nodal initial position  
            node_positions[i] = [
                mesh.NodeById(node).X, 
                mesh.NodeById(node).Y, 
                mesh.NodeById(node).Z
            ] 
    
    # Find node centroid and remove duplicate elements from list 
    centroid = avg_disp(node_positions) 
    elem_ids=list(set(elem_ids)) 

    # Count the number of nodes that the result should be summed over 
    # This will be n_nodes_per_element * n_elements 
    n_forces = get_n_elemnodal_forces(mesh, elem_ids, nodes)  

    # Retrieve the nodal forces and positions from the results file, then 
    # compute the moment in the specified local coordinate system 
    node_forces, node_positions = get_elemnodal_data(
        mesh, result_enfo, result_locdef, n_forces, elem_ids, nodes, is_ld_on, unit_scale
    )
    local_moment, r_max = calculate_moment(
        node_forces, node_positions, centroid, rotation_matrix
    )

    # Sign inverse for external force reaction 
    # (return the reaction force, which is opposite the internal force)
    local_moment = [-m for m in local_moment] 

    return local_moment, r_max


def process_section(analysis, body, local_csys, is_ld_on, unit_scale, model_scale):
    """ Compute the moment carried by a solid cross-section

    Args
    ---
    analysis: ANSYS ACT Analysis object
    body: ANSYS ACT GeoData Body object
    local_csys: ANSYS ACT Coordinate System object
        local coordinate system definition  
    is_ld_on : Bool
        True if NLGEOM is set to ON, False otherwise
    unit_scale: float 
        scale factor for unit system 
    model_scale: float 
        scale factor for model unit system 

    Returns
    ---
    [Mx, My, Mz] in global coordinate system

    """

    reader = analysis.GetResultsData()
    mesh = analysis.MeshData
    body_elements = mesh.MeshRegionById(body.Ids[0]).Elements
    n_body_elements = len(body_elements)
    result_locdef = reader.GetResult("LOC_DEF") 
    result_enfo = reader.GetResult("ENFO")
    rotation_matrix = rotation_matrix_from_local_csys(local_csys)
    origin = local_csys.Origin

    Min=[100000 for _ in range(n_body_elements)] 
    Max=[0 for _ in range(n_body_elements)] 

    # Elements in section and nodes on positive side of section (+Z)
    section_elements = [] 
    positive_nodes = [] 

    # Identify section elements containing nodes in +Z and -Z 
    for e, element in enumerate(body_elements): 
        for node in element.Nodes:
            r = [
                model_scale*node.X - origin[0], 
                model_scale*node.Y - origin[1], 
                model_scale*node.Z - origin[2]
            ]
            loc = transform(rotation_matrix, r)
            z_loc = loc[2]
            if z_loc < Min[e]: 
                Min[e] = z_loc
            if z_loc > Max[e]: 
                Max[e] = z_loc
            
        if Min[e] < 0 and Max[e] > 0: 
            section_elements.append(element)
            
    # Identify positive nodes from section elements
    for element in section_elements: 
        for node in element.Nodes:
            r = [
                model_scale*node.X - origin[0],
                model_scale*node.Y - origin[1],
                model_scale*node.Z - origin[2]]
            loc = transform(rotation_matrix, r)
            if loc[2] > 0:
                positive_nodes.append(node.Id)
                
        # Only include unique nodes 
        positive_nodes = list(set(positive_nodes))

    elem_ids = [elem.Id for elem in section_elements]
            
    n_forces = get_n_elemnodal_forces(mesh, elem_ids, positive_nodes)
    node_forces, node_positions = get_elemnodal_data(
        mesh, result_enfo, result_locdef, n_forces, elem_ids, positive_nodes, is_ld_on, unit_scale
    )
    centroid = avg_disp(node_positions)  
    local_moment, r_max = calculate_moment(
        node_forces, node_positions, centroid, rotation_matrix
    )
    
    # Sign inverse for external force reaction 
    # (return the reaction force, which is opposite the internal force)
    # TODO: should this be the same or opposite sign for section?
    local_moment = [m for m in local_moment] 

    return local_moment, r_max
    

def LDMProbe(result, stepInfo, collector):
    """ Main function, called by the tree object
    """

    # Initialize and collect objects and data
    analysis = result.Analysis
    reader = analysis.GetResultsData()
    reader.CurrentResultSet = stepInfo.Set
    solved_steps = reader.ListTimeFreq
    if stepInfo.Time not in solved_steps:
        return

    is_ld_on = analysis.AnalysisSettings.PropertyByName('UseLargeDeformation').StringValue == "On"
    unit_scale, model_scale = get_scale(analysis, is_ld_on)

    mode = result.Properties["Mode"].Value # "Interface" or "Section"
    local_csys = result.Properties["Orientation"].Value
    nodes = collector.Ids

    if mode == "Interface":
        local_moment, r_max = process_interface(analysis, nodes, local_csys, is_ld_on, unit_scale)

    elif mode == "Section":
        body = result.Properties["Geometry"].Value
        local_moment, r_max = process_section(analysis, body, local_csys, is_ld_on, unit_scale, model_scale)
        
    else:
        # TODO: print error message to user 
        # This should not be a condition that can happen...
        assert False 

    # Assign just the first node the result for display to the user
    collector.SetValues(nodes[0], [local_moment[0], local_moment[1], local_moment[2]])

    # Set initial graphics for triad and vector
    if result.DisplayTime.Value == 0: 
        # If Display Time is "Last", find numerical value. 
        check = max(solved_steps)
    else:
        check = result.DisplayTime.Value
        
    if stepInfo.Set == check: # Plot CS and M vector for selected time step display
        ExtAPI.Graphics.Scene.Clear()
        plot_csys_and_vector(result, mode, local_csys, local_moment, r_max, model_scale)

    reader.Dispose()


def plot_csys_and_vector(result, mode, local_csys, local_moment, r_max, model_scale):

    origin = local_csys.Origin

    # Store hidden result properties (see XML) for storing direction components
    result.Properties["Results/Mx"].Value = local_moment[0] 
    result.Properties["Results/My"].Value = local_moment[1]
    result.Properties["Results/Mz"].Value = local_moment[2]
    result.Properties["Results/Rmax"].Value = r_max

    if mode == "Interface":
        faceIds=result.Properties["Geometry"].Value
        centroids = []
        for id in faceIds:
            face = ExtAPI.DataModel.GeoData.GeoEntityById(id)
            centroids.append(face.Centroid)
        center = avg_disp(centroids)

    elif mode == "Section":
        center = [
            origin[0] / model_scale, 
            origin[1] / model_scale, 
            origin[2] / model_scale
        ]

    else:
        # TODO: print error message 
        assert False

    # Create triad graphics and set size 
    cs = ExtAPI.Graphics.Scene.Factory3D.CreateTriad(r_max) 

     # Convert 3x3 CS into STM 4x4 and define triad origin and orientation
    rotation_matrix = create_stm(local_csys.XAxis, local_csys.YAxis, local_csys.ZAxis)
    cs.Transformation3D.Set(rotation_matrix) 
    cs.Transformation3D.Translate(center[0],center[1],center[2]) 

    # Create vector graphics and set size 
    Vector = ExtAPI.Graphics.Scene.Factory3D.CreateArrow(2*r_max) 
    Vector.Color = 0x800080 # Purple
    x = result.Properties["Results/Mx"].Value
    y = result.Properties["Results/My"].Value
    z = result.Properties["Results/Mz"].Value

    # Set moment vector origin and orientation
    r = sqrt(y*y+z*z)
    Vector.Transformation3D.Set(rotation_matrix)
    # Vector defaults to Z axis, so rotate into result direction
    Vector.Transformation3D.Rotate(
        ExtAPI.Graphics.CreateVector3D(
            local_csys.YAxis[0], local_csys.YAxis[1], local_csys.YAxis[2]
        ), 
        atan2(x, r)
    ) 
    Vector.Transformation3D.Rotate(
        ExtAPI.Graphics.CreateVector3D(
            local_csys.XAxis[0], 
            local_csys.XAxis[1],
            local_csys.XAxis[2]
        ), 
        atan2(z,y) - pi/2.0
    )
    Vector.Transformation3D.Translate(center[0], center[1], center[2]) 


def ShowCS(result): 
    """ Display graphics when result object is selected in the tree 
    
    Does the following: 
    - Shows and scales vectors 
    - Converts results using proper model units 
    """
    
    if result.State == "solved":
        
        # Settings independent of tool mode
        ExtAPI.Graphics.ViewOptions.ShowResultVectors = True
        ExtAPI.Graphics.ViewOptions.VectorDisplay.DisplayType = VectorDisplayType.Sphere
        ExtAPI.Graphics.ViewOptions.VectorDisplay.LengthMultiplier = 0.1
        modunits = ExtAPI.DataModel.GeoData.Unit 
        scale = units.ConvertToUserUnit(ExtAPI, 1, modunits, "Length")
        mode = result.Properties["Mode"].Value
        CS = result.Properties["Orientation"].Value
        r_max = result.Properties["Results/Rmax"].Value
        
        if mode == "Interface":
            faceIds = result.Properties["Geometry"].Value
            centroids = []
            for id in faceIds:
                face = ExtAPI.DataModel.GeoData.GeoEntityById(id)
                centroids.append(face.Centroid)
            center = avg_disp(centroids)
            
        if mode == "Section":
            center = (CS.Origin[0]/scale, CS.Origin[1]/scale, CS.Origin[2]/scale)
            
        size = r_max
        cs = ExtAPI.Graphics.Scene.Factory3D.CreateTriad(size)
        M = create_stm(CS.XAxis, CS.YAxis, CS.ZAxis)
        cs.Transformation3D.Set(M)
        cs.Transformation3D.Translate(center[0],center[1],center[2])
        Vector = ExtAPI.Graphics.Scene.Factory3D.CreateArrow(2*size)
        Vector.Color = 0x800080 
        x = result.Properties["Results/Mx"].Value
        y = result.Properties["Results/My"].Value
        z = result.Properties["Results/Mz"].Value
        r = sqrt(y*y + z*z)
        Vector.Transformation3D.Set(M)
        Vector.Transformation3D.Rotate(ExtAPI.Graphics.CreateVector3D(CS.YAxis[0],CS.YAxis[1],CS.YAxis[2]), atan2(x,r))
        Vector.Transformation3D.Rotate(ExtAPI.Graphics.CreateVector3D(CS.XAxis[0],CS.XAxis[1],CS.XAxis[2]), atan2(z,y) - pi/2.0)
        Vector.Transformation3D.Translate(center[0],center[1],center[2])


def HideCS(result):
    """ Clear graphics object when tool is de-selected
    """ 
    
    ExtAPI.Graphics.Scene.Clear()


def geoCheck(result, prop):
    """ Checks to make sure the geometry scoping is correct
    
    Args
    ---
    result: 
    prop: property of the Geometry/callback
    
    Returns
    ---
        True, if mode is "Interface" and geometry is a face
        True, if mode is "Section" and geometry is a body
        False, otherwise
    """
    
    is_mode_interface = result.Properties["Mode"].Value == "Interface"
    is_mode_section = result.Properties["Mode"].Value == "Section"
    is_geom_face = ExtAPI.DataModel.GeoData.GeoEntityById(prop.Value.Ids[0]).Type == GeoCellTypeEnum.GeoFace
    is_geom_body = ExtAPI.DataModel.GeoData.GeoEntityById(prop.Value.Ids[0]).Type == GeoCellTypeEnum.GeoBody
        
    if prop.Value is not None: 
        if is_mode_interface and is_geom_face:
            return True 
        
        elif is_mode_section and is_geom_body:
            return True 
            
        else: 
            return False
    
    else:
        return False 