""" Large Displacement Moment Probe for ANSYS Mechanical
    (c) 2026 Yair Preiss

TODO: Short explanation of what the tool does, how it works, and how to install it.
"""

import units
import math
from System.Collections.Generic import List
from System import Array


#-------------------Helper Functions--------------------


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
    x, y, z: the 3-length rows of the spatial transformation matrix corresponding to the 
        direction of each axis in the global system
    
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


def get_unit_scale(reader, is_ld_on):
    model_units = ExtAPI.DataModel.GeoData.Unit 
    result_locdef = reader.GetResult("LOC_DEF")
    solve_units = result_locdef.GetComponentInfo('X').Unit
    model_scale = units.ConvertToUserUnit(ExtAPI, 1, model_units, "Length")

    if is_ld_on:
        unit_scale = units.ConvertUnit(1, solve_units, model_units,"Length")

    return unit_scale


def get_n_elemnodal_forces(mesh, elem_ids, node_ids):
    count = 0 

    for elem_id in elem_ids:
        element = mesh.ElementById(elem_id)
        for node_id in node_ids:
            if element.NodeIds.IndexOf(node_id) >= 0: 
                # Nodes attached to the element have a positive index
                count += 1
    
    return count


def get_elemnodal_data(mesh, result_enfo, result_locdef, n_forces, elem_ids, node_ids, is_ld_on, unit_scale):
    node_forces = [[] for _ in range(n_forces)]  
    node_positions = [None for _ in range(n_forces)]

    count = 0
    for Id in elem_ids:
        element = mesh.ElementById(Id)
        # elementnodal force reactions are reported in a row-major linear array, 
        # with each [Fx1, Fy1, Fz1, Fx2, Fy2, Fz3, ... Fzn], not as individual vectors
        elem_force = result_enfo.GetElementValues(Id) 
        index = List[object]()

        # TODO: this is repeated work, we could just associate node ids with the element id
        for node_id in node_ids:
            if element.NodeIds.IndexOf(node_id) >= 0:
                index.Add(element.NodeIds.IndexOf(node_id))

        for i in range(len(index)):
            # Assign elementnodal force reaction into node-specific vector
            node_forces[count] = [0.0 for _ in range(3)]
            node_forces[count][0] = elem_force[3*index[i]]
            node_forces[count][1] = elem_force[3*index[i]+1]
            node_forces[count][2] = elem_force[3*index[i]+2]
            if is_ld_on:
                # Assign elementnodal position vectors to pair with respective force reactions
                node_positions[count] = [x*unit_scale for x in result_locdef.GetNodeValues(element.NodeIds[index[i]])] 
            else:
                node_positions[count] = [
                    mesh.NodeById(element.NodeIds[index[i]]).X, 
                    mesh.NodeById(element.NodeIds[index[i]]).Y, 
                    mesh.NodeById(element.NodeIds[index[i]]).Z
                ]
            count += 1

    return node_forces, node_positions


def process_interface(analysis, nodes, is_ld_on, unit_scale, rotation_matrix):
    """ If mode is interface, process the reaction moment
    """

    # Unpack data related to the model 
    reader = analysis.GetResultsData()
    result_locdef = reader.GetResult("LOC_DEF") 
    result_enfo = reader.GetResult("ENFO")
    mesh = analysis.MeshData 
    
    n_nodes = len(nodes)

    # Node positions; this will become a list of 3-length lists
    node_positions=[None for _ in range(n_nodes)] 

    # Number of element ids is not known a priori 
    elem_ids = []

    # Populate node positions array
    for i, node in enumerate(nodes):
        # Get all elements associated with nodes
        elem_ids = elem_ids + [int(x) for x in mesh.NodeById(node).ConnectedElementIds] 
        if is_ld_on:
            # Get nodal positions at load step
            node_positions[i] = [x * unit_scale for x in result_locdef.GetNodeValues(node)] 
        else:
            # Get nodal initial position  
            node_positions[i] = [mesh.NodeById(node).X, mesh.NodeById(node).Y, mesh.NodeById(node).Z] 
    
    # Find node centroid and remove duplicate elements from list 
    centroid = avg_disp(node_positions) 
    elem_ids=list(set(elem_ids)) 

    # Count the number of nodes that the result should be summed over 
    # This will be n_nodes_per_element * n_elements 
    n_forces = get_n_elemnodal_forces(mesh, elem_ids, nodes)  

    node_forces, node_positions = get_elemnodal_data(mesh, result_enfo, result_locdef, n_forces, elem_ids, nodes, is_ld_on, unit_scale)

    # Compute relative position after deflection (local r vector)
    # Also keep track of largest value for sizing the display vector
    r = [None for _ in range(n_forces)] 
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
        global_moment = vsum(global_moment,cross(r[f],force)) 

    local_moment = transform(rotation_matrix,global_moment) # Transform M into selected CS orientation 

    # Sign inverse for external force Reaction 
    # (return the reaction force, which is opposite the internal force)
    local_moment = [-m for m in local_moment] 


#-------------------Main Fuction--------------------


def LDMProbe(result,stepInfo,collector):

    #Initialize and collect objects and data
    analysis = result.Analysis
    reader = analysis.GetResultsData()
    solved_steps=reader.ListTimeFreq
    if stepInfo.Time not in solved_steps:
        return
    mesh = analysis.MeshData

    is_ld_on = analysis.AnalysisSettings.PropertyByName('UseLargeDeformation').StringValue == "On"
    unit_scale = get_unit_scale(reader, is_ld_on)

    mode = result.Properties["Mode"].Value # "Interface" or "Section"
    CS = result.Properties["Orientation"].Value
    reader.CurrentResultSet = stepInfo.Set
    nodes = collector.Ids
    rotation_matrix = [CS.XAxis, CS.YAxis, CS.ZAxis]

    #Begin processing
    if mode == "Interface":
        process_interface(analysis, nodes, is_ld_on, unit_scale, rotation_matrix)
    if mode == "Section":
        body = result.Properties["Geometry"].Value
        totelement = mesh.MeshRegionById(body.Ids[0]).Elements
        Min=[100000 for l in range(len(totelement))] # Initialize Min check value
        Max=[0 for l in range(len(totelement))] # Initialize Max check value
        secelement = [] # Elements in section
        posnode = [] # Nodes on positive side of section (+Z)
        for e,element in enumerate(totelement): # Identify section elements containing nodes in +Z and -Z 
            for node in element.Nodes:
                o = [mod_scale*node.X-CS.Origin[0],mod_scale*node.Y-CS.Origin[1],mod_scale*node.Z-CS.Origin[2]]
                loc = transform(rot,o)
                z_loc = loc[2]
                if z_loc<Min[e]: Min[e]=z_loc
                if z_loc>Max[e]: Max[e]=z_loc
            if Min[e]<0 and Max[e]>0: 
                secelement.Add(element)
        for element in secelement: # Identify positive nodes from section elements
            for node in element.Nodes:
                o = [mod_scale*node.X-CS.Origin[0],mod_scale*node.Y-CS.Origin[1],mod_scale*node.Z-CS.Origin[2]]
                loc = transform(rot,o)
                if loc[2]>0:
                    posnode.Add(node.Id)
            posnode = list(set(posnode))
        n = len(posnode)
        n_pos=[None]*n 
        for i, node in enumerate(posnode):
            if LD=="On":
                n_pos[i] = [x*sol_scale for x in LOCDEF.GetNodeValues(node)] # Get nodal positions at load step
            if LD=="Off":
                n_pos[i] = [mesh.NodeById(node).X, mesh.NodeById(node).Y, mesh.NodeById(node).Z] # Get nodal initial position
        CG = avg_disp(n_pos) # Determine nodal centroid            
        F_count = 0
        for element in secelement:
            for node in posnode:
                if element.NodeIds.IndexOf(node)>=0:
                    F_count=F_count+1
        F_nodes = [None]*F_count
        n_pos = [None]*F_count
        count = 0
        for element in secelement:
            F_element = F.GetElementValues(element.Id) # elementnodal force reactions are reported in one long list with each [Fx1, Fy1, Fz1, Fx2, Fy2, Fz3, ... Fzn], not as individual vectors
            index = List[object]()
            for node in posnode:
                if element.NodeIds.IndexOf(node)>=0:
                    index.Add(element.NodeIds.IndexOf(node))
            for i in range(len(index)):
                F_nodes[count]=[None]*3 # Assign elementnodal force reaction into node-specific vector
                F_nodes[count][0] = F_element[3*index[i]]
                F_nodes[count][1] = F_element[3*index[i]+1]
                F_nodes[count][2] = F_element[3*index[i]+2]
                if LD=="On":
                    n_pos[count] = [x*sol_scale for x in LOCDEF.GetNodeValues(element.NodeIds[index[i]])] # Assign elementnodal position vectors to pair with respective force reactions
                if LD=="Off":
                    n_pos[count] = [mesh.NodeById(element.NodeIds[index[i]]).X, mesh.NodeById(element.NodeIds[index[i]]).Y, mesh.NodeById(element.NodeIds[index[i]]).Z]                
                count = count+1  
        r = [None]*F_count 
        Global_M = [0.0,0.0,0.0]
        Local_M = [0.0,0.0,0.0]
        r_max = 0
        for f, force in enumerate(F_nodes):
            r[f] = [(n_pos[f][0]-CG[0]),(n_pos[f][1]-CG[1]),(n_pos[f][2]-CG[2])] # Determine loacal r vector
            if mag(r[f]) > r_max:
                r_max = mag(r[f])
            Global_M = vsum(Global_M,cross(r[f],force)) # Sum all M=rxF products in Global orientation (default solver reporting)
        Local_M = transform(rot,Global_M) # Transform M into selected CS orientation 
        Local_M = [m for m in Local_M] # Sign inverse for external force Reaction.       

    collector.SetValues(nodes[0],[Local_M[0],Local_M[1],Local_M[2]])

    #Set initial graphics for triad and vector
    if result.DisplayTime.Value == 0: # If Display Time is "Last", find numerical value. 
        check = max(solved_steps)
    else:
        check = result.DisplayTime.Value
    if stepInfo.Set == check: # Plot CS and M vector for selected time step display
        ExtAPI.Graphics.Scene.Clear()
        result.Properties["Results/Mx"].Value = Local_M[0] # Store hidden result properties (see XML) for storing direction components
        result.Properties["Results/My"].Value = Local_M[1]
        result.Properties["Results/Mz"].Value = Local_M[2]
        result.Properties["Results/Rmax"].Value = r_max
        if mode == "Interface":
            faceIds=result.Properties["Geometry"].Value
            centroids=[]
            for id in faceIds:
                face = ExtAPI.DataModel.GeoData.GeoEntityById(id)
                centroids.append(face.Centroid)
            center = avg_disp(centroids)
        if mode == "Section":
            center = (CS.Origin[0]/mod_scale,CS.Origin[1]/mod_scale,CS.Origin[2]/mod_scale)
        size = r_max
        cs = ExtAPI.Graphics.Scene.Factory3D.CreateTriad(size) # Create triad graphics and set size 
        M=create_stm(CS.XAxis,CS.YAxis,CS.ZAxis) # Convert 3x3 CS into STM 4x4
        cs.Transformation3D.Set(M) # Define triad orientation
        cs.Transformation3D.Translate(center[0],center[1],center[2]) # Define tirad origin
        Vector = ExtAPI.Graphics.Scene.Factory3D.CreateArrow(2*size) # Create vector graphics and set size 
        Vector.Color = 0x800080 # Purple
        x = result.Properties["Results/Mx"].Value
        y = result.Properties["Results/My"].Value
        z = result.Properties["Results/Mz"].Value
        r = sqrt(y*y+z*z)
        Vector.Transformation3D.Set(M)  # Define vector orientation
        Vector.Transformation3D.Rotate(ExtAPI.Graphics.CreateVector3D(CS.YAxis[0],CS.YAxis[1],CS.YAxis[2]), atan2(x,r)) # Vector defaults to Z axis, so rotate into result direction
        Vector.Transformation3D.Rotate(ExtAPI.Graphics.CreateVector3D(CS.XAxis[0],CS.XAxis[1],CS.XAxis[2]), atan2(z,y) - pi/2.0)
        Vector.Transformation3D.Translate(center[0],center[1],center[2]) # Define vector origin
    reader.Dispose()


#-------------------Graphics Display/Hide--------------------


def ShowCS(result): # Graphics to display when result object is selected in Tree
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


#-------------------Scoping verification--------------------


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
        if is_mode_interface:
            if is_geom_face:
                return True 
            else: 
                return False 
        
        elif is_mode_section:
            if is_geom_body:
                return True 
            else:
                return False 
        
        else: 
            return False 
    
    else:
        return False 