# pocket_capping

### Description

Using a .mol2 protein-complex, this script selects the pocket residues and cap incomplete residues with hydrogen atoms.

First, the script finds and select residues within a cutoff radius that is defined by the user. The pocket selection will include water molecules and other residues within the cutoff.

Then, using graph theory and a residue graph dictionary as a reference, the script finds which residues will not be entirely connected. With this information and the protein-complex as a reference, the script caps those residues with hydrogen in their corresponding position and at an equilibrium distance that can be defined by the user.

Finally, the script will output a .mol2 file of the pocket capped with hydrogens. For the moment the output file does not have bond information.

For the moment the script can only cap the residues with hydrogen, but it will be updated to cap residues with ACE and NME or with user-defined entities.

| ![](https://github.com/user-attachments/assets/eb2a4d35-86f1-4f53-8fba-4b577259fdb1) | 
|:--:| 
| *Pocket capped with hygrogen atoms. Yellow: the ligand; blue: connected residues; red: disconnected residues; green: hydrogens added; white: the protein.* |
