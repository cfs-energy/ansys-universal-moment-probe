# Ansys NLGEOM Moment Probe Extension

[Full Manual](Docs/Manual.md)  

An add-on for Ansys Mechanical for a Moment Probe result object that self-corrects for Large Deflection (NLGEOM = ON) effects.

## Installation
This extension can be installed in two ways: binary installer, or directory association.
### Binary Installation
Do one of the following:
1. Workbench menu bar:  
1.1. Extensions > Install Extension...  
1.2. Select the binary (*.wbex) file.  
1.3. Extensions > Manage Extensions...  
1.4. Click the check box next to the extension you want to activate.  
2. ACT Start Page:  
2.1. Below the Workbench menu bar, click on "ACT Start Page".  
2.1. Go to "Manage Extensions" (middle of ACT Start Page window).  
2.3. In the upper right corner, click the "+" button.  
2.4. Select the binary file.  
2.5. The new extension should be added as a grey block into the collection of extensions. Simply click on it to load it (will unload after you close Workbench), or click on the upside-down triangle on the bottom right corner of the extension block and select "Load as default" to have it permanently loaded.

### Directory Association
1. Download the contents of this repo into a dedicated folder
2. Below the Workbench menu bar, click on "ACT Start Page".  
3. Go to "Manage Extensions" (middle of ACT Start Page window).
4. In the upper right corner, click the Gear button.  
5. Click the "+ Add Folder" button.
6. Select the directory where the LMDProbe.xml file is located. 
7. The new extension should be added as a grey block into the collection of extensions. Simply click on it to load it (will unload after you close Workbench), or click on the upside-down triangle on the bottom right corner of the extension block and select "Load as default" to have it permanently loaded.

**NOTE**: When working with directory association, you can modify the code directly and simply click the "Reload" button (Mechanical > Automation tab > ACT Development block) to update the loaded extension.


## Contributing
Contributions in the way of bugfixes and new features are welcome!  

For requests or questions, please contact [Yair Preiss](mailto:ypreiss@cfs.energy).
