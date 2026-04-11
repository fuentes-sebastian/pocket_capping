# pocket_capping
From a .mol2 protein-complex, this script selects the pocket residues and cap incomplete residues.

First, the script finds and select residues within a cutoff radius that is defined by the user. Then, using graph theory, the script finds which residues will not be completely connected. Using the protein-complex as a reference, the script caps those residues with hydrogen in the "correct" position and at an "adequate" equilibrium distance that can also be defined by the user.

For the moment the script can only cap the residues with hydrogen, but it will be updated to cap residues with ACE and NME or with user-defined entities.
<figure>
<img width="546" height="576" alt="Capping_example" src="https://github.com/user-attachments/assets/eb2a4d35-86f1-4f53-8fba-4b577259fdb1" 
  <figcaption align="center"> <b> Figure 1: </b> Example of the capping. The ligand is coloured in yellow, blue residues are bonded residues while red residues are "incomplete" to-be-capped residues. Green spheres show the hydrogens added to the uncapped residues and white lines represent the original protein.> 
<\figure>
