# pocket_capping

### Description

Using a .mol2 protein-complex, this script selects the pocket residues and cap incomplete residues with hydrogen atoms.

First, the script finds and select residues within a cutoff radius that is defined by the user. 

Then, using graph theory, the script finds which residues will not be entirely connected. Using the protein-complex as a reference, the script caps those residues with hydrogen in the corresponding position and at an equilibrium distance that can be defined by the user.

For the moment the script can only cap the residues with hydrogen, but it will be updated to cap residues with ACE and NME or with user-defined entities.

![Example of the capping. The ligand is coloured in yellow, blue residues are bonded residues while red residues are "incomplete" to-be-capped residues. Green spheres show the hydrogens added to the uncapped residues and white lines represent the original protein.](https://github.com/user-attachments/assets/eb2a4d35-86f1-4f53-8fba-4b577259fdb1)
