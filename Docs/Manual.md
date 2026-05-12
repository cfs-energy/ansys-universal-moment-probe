# Extension Manual

## Purpose
When Large Deflection is on, the standard ANSYS Moment Probe will calculate the moment about the undeformed position of the scoped geometry, with the nodal forces at the displaced nodal positions. This means that the resultant moments will have fictitious moment arms and the values incorrect. This installable adds a Tree Object result into the Solution right-click menu, designed to replicate the functionality of the built-in Moment Probe, but with compensation for nodal displacements when Large Deflection is set to “On”.  

## Instrucitons for Use 

### Result Value Data

### A Note on Graphics

## Core Process Outline
### Interface Mode
1. Nodes collected from scoped Face.  
2. Location of each node is collected and used to determine centroid:  
    2.1 Large Deflection = “On”: uses “LOC_DEF”.  
    2.2 Large Deflection = “Off”, uses initial Node position.  
3. All Elements associated with Nodes are collected.   
4. Element-Nodal reaction loads (“ENFO”) for each Element are collected (F).  
5. For each Node (from nodes collected in [1.]) within each Element, the position vector (r) is determined with respect to the centroid of the Nodes:  
    5.1  Large Deflection = “On”: uses “LOC_DEF”.  
    5.2 Large Deflection = “Off”, uses initial Node position.  
6. For each Node* within each Element, a moment (M) is calculated about the displaced centroid by M=rxF. 
7. The total moment is then determined by vector sum of all nodal moments, and rotated to match selected Coordinate System.

*A Node is considered once for every Element (and its internal force reaction) it is associated with, meaning that it might be invoked several times in this step.   

**NOTE**: Orientation coordinate system is not used as the point of moment summation, but only for orientation. This means the selected CS does not have to be located at the target geometry, and that the same CS can be used for any number of probes as long as the desired orientation is common.

**NOTE**: The output data table will include a Maximum, Minimum, and Total columns for each component. These values will all the identical and the Total column should not be confused with a vector sum.  

### Section Mode
This mode follows a similar process to the Interface Mode, except that the nodes used are filtered down in the following manner:
1. Obtain all elements in the scoped body.
2. Define section plane based on scoped CS Z-axis normal.
3. Filter down only elements in section plane.
4. Filter down only element nodes nodes on one side of the section plane (positive Z-axis).

**NOTE**: Orientation coordinate system is not used as the point of moment summation, but only for section location and result orientation.

**NOTE**: The output data table will include a Maximum, Minimum, and Total columns for each component. These values will all the identical and the Total column should not be confused with a vector sum. 

## Example Case

The example discussed here will not elaborate on the original Ansys moment probe (that does not correct for NLGEOM), but on the differences between this extention and the new Ansys moment probe (as of 2026R1) with the NLGEOM correction active. The model presented here is intented as an example of the error that can occur when using the Ansys probe, which is in no way limited to this specific condition. 

### Model Setup and Results

In this model we have a hollow tube fixed at one end and with a force aplied at the other (noted refernce to the CS shown in the image below). A coordinate system for probing moments is placed 10mm from the end of the tube. Two load cases are considered:  
1. Fz=1000N
2. Fz=Fx=1000N    

<p align="center" width="100%">
  <img src="img/pipe.png" alt="radius1" width="75%" align="center"/>
</p>

While for the first load case with only the axial load we will expect no moment in any axis, the second case requires some evaluation by hand calcs. The following two images show the displacement in the X-Axis (top) and Z-Axis at the location of the probe and at the face where the force is being applied.  

<p align="center" width="100%">
  <img src="img/disp.png" alt="radius1" width="75%" align="center"/>
</p>

Knowing that the initial moment arm for the Fx component was 10mm, and 0mm for the Fz, we can calculated the resultant moment for that section:
$$
M = F_x \times (r_z - (dz_{tip} - dz_{section})) - F_z \times (r_x - (dx_{tip} - dx_{section}))
$$ 
and simplifying for Fx=Fz, we get:
$$
M = 1000[N] \times (10 - (3.426125765 - 3.041292906) - (30.52874565 - 27.80307388))[mm] = 6,889[Nmm]
$$ 

Now we can compare these values against the Ansys moment probe and this extesion:
<div align="center">

| Load Case | Hand Calc | Ansys Probe | This Extension |
|-----------|-----------|-------------|----------------|
|1          |0 [Nmm]    |4,500 [Nmm]  |0.5 [Nmm]       |
|2          |6,889[Nmm] |2,004 [Nmm]  |6,980 [Nmm]     |

</div>

### Discussion 

The reason for the error in the Ansys probe is that it calculates the moment about the position (LOCDEF) of the nearest node to the scoped CS, not the CS itself, when NLGEOM=on (as shown below). In the first load step the purely axial load of 1,000 N is calculated at a node on the ID of the tube (r=4.5 mm), resulting in 4,500 Nmm of ficticious moment. This is also true for the second load case, though it is less obvious since the moment is not expected to be 0.    

<p align="center" width="100%">
  <img src="img/ansysprobe.png" alt="radius1" width="50%" align="center"/>
</p>

It should be noted than each time this particualr model is re-meshed and re-run, the total moment on the second step as measured by the Ansys probe reslted in a different value. This is because the "closet node" among all equidistant nodes on the ID, with respect to the center CS, is randomly assigned, and so the center of moment summation is not constant.

### Conclusion

More generally, it should be observed that the Ansys moment probe will always be incorrect unless there happens to be a mesh node at the exact position of the desired CS. Error will be proportional to the distance to the closest node and the magnitude of load perpendicular to that distance vector.  

**CRITICAL NOTE**: The Ansys moment probe output will not necessarily be higher than the true value, but it will be wrong. Using the Ansys moment probe does NOT guarantee conservatism.