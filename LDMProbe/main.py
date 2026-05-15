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
    

def create_create_stm(x, y, z):  
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
        
    matrix[3]=0.0
    matrix[7]=0.0
    matrix[11]=0
    matrix[15]=1.0
    
    return matrix        


def cross(a, b):     
    """ Return the cross product of two 3-length lists (vectors)
    """
    
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
         
    return c

#Coordinate transformation
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
    disp = [sum(coord) for coord in zip(*points)] / n
    return disp


#-------------------Main Fuction--------------------


def LDMProbe(result,stepInfo,collector):

    #Initialize and collect objects and data
    analysis = result.Analysis
    reader = analysis.GetResultsData()
    solved_steps=reader.ListTimeFreq
    if stepInfo.Time not in solved_steps:
        return
    mesh = analysis.MeshData
    modunits = ExtAPI.DataModel.GeoData.Unit
    LD=analysis.AnalysisSettings.PropertyByName('UseLargeDeformation').StringValue
    LOCDEF = reader.GetResult("LOC_DEF") 
    F = reader.GetResult("ENFO")
    mode = result.Properties["Mode"].Value
    CS = result.Properties["Orientation"].Value
    solve_unit = LOCDEF.GetComponentInfo('X').Unit
    mod_scale = units.ConvertToUserUnit(ExtAPI,1,modunits,"Length")
    if LD=="On": 
        sol_scale = units.ConvertUnit(1,solve_unit,modunits,"Length") # Set scale factor between Solver Units and Model Units
    reader.CurrentResultSet = stepInfo.Set
    nodes = collector.Ids
    rot = [CS.XAxis,CS.YAxis,CS.ZAxis]

    #Begin processing
    if mode == "Interface":
        n = len(nodes)
        n_pos=[None]*n 
        elementIds = []
        for i, node in enumerate(nodes):
            elementIds = elementIds+[int(x) for x in mesh.NodeById(node).ConnectedElementIds] # Get all elements associated with nodes
            if LD=="On":
                n_pos[i] = [x*sol_scale for x in LOCDEF.GetNodeValues(node)] # Get nodal positions at load step
            if LD=="Off":
                n_pos[i] = [mesh.NodeById(node).X, mesh.NodeById(node).Y, mesh.NodeById(node).Z] # Get nodal initial position  
        CG = avg_disp(n_pos) # Determine nodal centroid
        elementIds=list(set(elementIds)) # Remove duplicate elements from list 
        F_count = 0
        for Id in elementIds:
            element = mesh.ElementById(Id)
            for node in nodes:
                if element.NodeIds.IndexOf(node)>=0: # Index will be -1 for nodes NOT in element
                    F_count = F_count+1 # Count number of elementnodal forces needed to resolve moment (each node, once per element)
        F_nodes = [None]*F_count
        n_pos = [None]*F_count
        count = 0
        for Id in elementIds:
            element = mesh.ElementById(Id)
            F_element = F.GetElementValues(Id) # elementnodal force reactions are reported in one long list with each [Fx1, Fy1, Fz1, Fx2, Fy2, Fz3, ... Fzn], not as individual vectors
            index = List[object]()
            for node in nodes:
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
        r_max = 0
        for f, force in enumerate(F_nodes):
            r[f] = [(n_pos[f][0]-CG[0]),(n_pos[f][1]-CG[1]),(n_pos[f][2]-CG[2])] # Determine loacal r vector
            if mag(r[f]) > r_max:
                r_max = mag(r[f])
            Global_M = vsum(Global_M,cross(r[f],force)) # Sum all M=rxF products in Global orientation (default solver reporting)
        Local_M = transform(rot,Global_M) # Transform M into selected CS orientation 
        Local_M = [-m for m in Local_M] # Sign inverse for external force Reaction.
    
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
    """ Checks ??? to make sure the scoping is correct
    """
    
    is_mode_interface = result.Properties["Mode"].Value == "Interface"
    is_mode_section = result.Properties["Mode"].Value == "Section"
    is_geom_face = ExtAPI.DataModel.GeoData.GeoEntityById(prop.Value.Ids[0]).Type == GeoCellTypeEnum.GeoFace
    is_geom_body = ExtAPI.DataModel.GeoData.GeoEntityById(prop.Value.Ids[0]).Type == GeoCellTypeEnum.GeoFace
        
    if prop.Value != None: 
        if is_mode_interface:
            if is_geom_face:
                return True 
            else: 
                return False 
        
        elif is_mode_section:
            if is_geom_face:
                return True; 
            else:
                return False 
        
        else: 
            return False 
    
    else:
        return False 