import networkx as nx # for graph comparison
import ast # to convert txt file to a list of dictionaries 
import numpy as np

# NOTE: this script searchs for a pocket, 
# identify pocket residues that will be incomplete using graph theory and
# cap incomplete residues with H using the complete protein as reference
# finally, output the coordinates of the pocket with the ligand in .mol2 format

# ABOUT INPUT FILE: must be a .mol2 file 
# Input file MUST have all hidrogens
# should also work with .itp/.top, but those are not tested

# ___________________________ variables ___________________________ 
input_file='input_file.mol2' # must have all hydrogens and bonds, input file

# pocket finding
cutoff=15 # Angstroms
lig_name='UNK' # ligand name (in this case, will search residue that starts with LIG name, dont have to be exact match)
r_CH=1.09 # distance of H capping the C atom, in Angstroms
r_NH=1.06 # distance of H capping the N atom, in Angstroms

# output file
output_file='pocket_test.mol2'
mol_name='molecule' # name of molecule (for .mol2 output file)

# Weight dictionary, needed to parse mol2 files 
weight_dict={'H':1, 'C':12, 'N':14, 'O':16, 'S':32}


# ___________________________ file reading ___________________________ 
warnings=[]

# reading file dictionaries (must be a list of dictionaries)
with open('dictionaries/dict_ref.txt', 'r', encoding='utf-8') as f:
    ref_dict = ast.literal_eval(f.read())

# reading system (.mol2, .itp or .top) file
with open(f"{input_file}", 'r', encoding="utf-8") as f:
    lines = f.readlines()

# ___________________________ parsing of .itp or .top file ___________________________
def parse_atoms(line):
    if not line.startswith(';'):
        atom=line.split()
        atom_info = {
            'id': int(atom[0]),
            'element': atom[4][0],
            'res_id': int(atom[2]),
            'res_name': atom[3],
            'weight':round(float(atom[7])),
            'coords':np.asarray([atom[2], atom[3], atom[4]])}
        return(atom_info)

def bonds_prase(line):
    if not line.startswith(';'):
        bond=line.split()
        bond_info = {
            'a_i': int(bond[0]),
            'a_j': int(bond[1])}
        return(bond_info)

# ___________________________ parsing .mol2 files ___________________________

def parse_atoms_mol2(line):
    atom=line.split()
    element=atom[1][0]

    if element not in weight_dict.keys():
        weight=0
        warn=f'atom with id {int(atom[0])} has element {element} which is not in dictionay of weights, wheight = 0 applied'
        print(warn)
        warnings.append(warn)
    if element in weight_dict.keys():
        weight=weight_dict[element]

    atom_info = {
        'id': int(atom[0]),
        'atom_name':atom[1],
        'element': element,
        'res_id': int(atom[6]),
        'res_name': atom[7] if len(atom) > 6 else 'Unk',
        'weight':weight, 
        'coords':np.asarray([float(atom[2]), float(atom[3]), float(atom[4])])}
    return(atom_info)

def bonds_prase_mol2(line):
    bond=line.split()
    bond_info = {
        'id': int(bond[0]),
        'a_i': int(bond[1]),
        'a_j': int(bond[2]),
        'bond_type' : int(bond[3])}
    return(bond_info)

# ___________________________ .mol2 file writing ___________________________
def mol_atom_writer(a_id, a_name, a_coords, a_element, r_id, r_name):
    line=(str(a_id).ljust(9)+
    str(a_name).rjust(6)+
    str('%8.3f' % a_coords[0]).ljust(9)+
    str('%8.3f' % a_coords[1]).ljust(9)+
    str('%8.3f' % a_coords[2]).ljust(9)+
    str(a_element).ljust(7)+
    str(r_id).ljust(7)+
    str(r_name).ljust(10))
    
    return(line)

def write_mol2(output_file, mol_name, n_atoms, n_res, atom_lines):
    with open(output_file, 'w') as f:
    
            # molecule
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{mol_name}\n")
            # num_atoms, num_bonds, num_substructures (usually 1 or 0), num_features, num_sets
            f.write(f"{str(n_atoms).ljust(8)} {str(0).ljust(8)} {str(n_res).ljust(8)}\n")
            f.write("PROTEIN\n")  # Molecule type (can also be SMALL or PROTEIN)
            f.write("USER_CHARGES\n") # Charge type (can be NO_CHARGES, GASTEIGER, etc.)
            
            # atoms
            f.write("@<TRIPOS>ATOM\n")
            for l in atom_lines:
                f.write(l + '\n')


# ___________________________ determine non-overlapping (unique) matches ___________________________
def find_uniques_matches(matches_list, total_nodes, ignore_total_nodes=True):
    # if matching doesn't contain all atoms of the residue, change ignore_total_nodes=True

    # Create a compatibility graph to determine overlapping
    compatibility_graph = nx.Graph()

    # Add nodes to compatibility graph
    for i in range(len(matches_list)):
        compatibility_graph.add_node(i)

    # Connect two matches if they don't overlap
    for i in range(len(matches_list)):
        nodes_i = set(matches_list[i].keys()) # for each match, a set of its atoms ids

        for j in range(i + 1, len(matches_list)):
            nodes_j = set(matches_list[j].keys()) # for each other match, other set of its atoms ids

            # If they don't share atoms. Connect them.
            # Makes a connection between each match that with different atoms
            if nodes_i.isdisjoint(nodes_j):
                compatibility_graph.add_edge(i, j)

    # Find cliques (overlapping)
    # Clique is group of nodes interconnected to everyone.
    # This means that all matched in cliques have diffetent atoms (thus, no overlapping)
    all_cliques = list(nx.find_cliques(compatibility_graph))

    # Convert indices back to matching data
    solutions = []
    for clique in all_cliques:
        solution_set = [matches_list[i] for i in clique]
        solutions.append(solution_set)

    # filter solutions so n_nodes = n_total_nodes
    if ignore_total_nodes==False:
        final_solutions=[]
        for sol in solutions:
            n_nodes = sum(len(m) for m in sol) # count nodes in solutions
            if n_nodes==total_nodes:
                final_solutions.append(sol)
        solutions=final_solutions
    return solutions

# ___________________________ graph creation ___________________________
def create_graph(node_list, name_list, edge_list, weight_list):
    graph_ref = nx.Graph()
    dict_ref=[]
    graph_ref.add_edges_from(edge_list)
    for n, nm, w in zip(node_list, name_list,weight_list):
        dict_ref.append({'id': n, 'name':nm})
        graph_ref.add_edge(n,-w) # adds atom mass as an edge between atom and its mass (we use neg value to avoid atom id matching with weight node)
    return(graph_ref, dict_ref)


# ___________________________ file parsing ___________________________ 
atom_data=[]
bond_data=[]
section=None

# parsing if .top or . itp
if input_file.endswith('.top') or input_file.endswith('.itp'):
    for line in lines:
        line = line.strip() # eliminate space between lines
        if line: # ignores empty lines
            if line.startswith("[ "):
                section = line
                continue # this ignores section line
            
            # Atoms data prasing
            if section == "[ atoms ]":
                a_info=parse_atoms(line)
                if a_info:
                    atom_data.append(a_info)
            
            if section == "[ bonds ]":
                b_info=bonds_prase(line)
                if b_info:
                    bond_data.append(b_info)

# parsing if .mol2
if input_file.endswith('.mol2'): 
    for line in lines:
        line = line.strip() # eliminate space between lines
        if line.startswith("@<TRIPOS>"): # this ignores section line
            section = line
            continue

        # Atoms data prasing
        if section == "@<TRIPOS>ATOM":
            a_info=parse_atoms_mol2(line)
            atom_data.append(a_info)

        if section == "@<TRIPOS>BOND":
            bond=line.split()
            b_info=bonds_prase_mol2(line)
            bond_data.append(b_info)
else:
    warn=f'File format for {input_file} not recognized'
    warnings.append(warn)
    print(warn)


#  ___________________________ Finding pocket  ___________________________ 
# extract ligand Data
lig_dict=[]

lig_coords=[]
lig_atom_ids=[]
lig_atom_names=[]
lig_elements=[]
for a in atom_data:
    # if a['res_name'] == lig_name: # use this for exact ligand name
    if a['res_name'].startswith(f'{lig_name}'):
        lig_coords.append(a['coords'])
        res_id=a['res_id']
        lig_name=a['res_name']
        lig_atom_ids.append(a['id'])
        lig_atom_names.append(a['atom_name'])
        lig_elements.append(a['element'])

lig_dict.append({'res_id':res_id,
                 'res_name':lig_name,
                 'atom_ids':lig_atom_ids,
                 'atom_names':lig_atom_names,
                 'elements':lig_elements,
                 'coords':lig_coords
                 })            

lig_coords=np.asarray(lig_coords) # ligand coordinates

# selection of residues within cutoff (pocket residues)
pocket_res_ids=[]
for a in atom_data:
    if a['res_name'].startswith(f'{lig_name}') == False: # to avoit comparison with ligand atoms
        for l_coord in lig_coords:
            a_coord =a['coords']
            dist = np.linalg.norm(a_coord-l_coord)
            if dist <= cutoff:
                if a['res_id'] not in pocket_res_ids: # res names without repetition
                    pocket_res_ids.append(a['res_id'])


#  ___________________________ pocket Graph ___________________________ 
n=0
dict_res=[]
matches_in_total=[]
pocket_ids=[] # complete list of pocket atoms ids (to identify bonds within pocket)
pocket_dict=[]

for res_id in pocket_res_ids:
    graph_res = nx.Graph()
    graph_ref = nx.Graph() # This, I think, avoids error if first residue is not in library

    res_atom_ids=[] # nodes created for each residue
    res_bonds=[] # list of bonds inside each residue
    res_dict=[] # dictionaries of atoms inside each residue
    res_weights=[]
    elements_list=[]
    coords_list=[]
    atom_names_list=[]
    
    for a in atom_data:

        # extraction of pocket data
        if a['res_id']==res_id:
            res_name_test=a['res_name']
            res_dict.append(a)
            res_atom_ids.append(a['id'])
            res_weights.append(a['weight'])
            elements_list.append(a['element'])
            atom_names_list.append(a['atom_name'])
            coords_list.append(a['coords'])
            pocket_ids.append(a['id'])
            n=a['res_id']
        
    # dictionary of pocket residues
    pocket_dict.append({'res_id':n,
                        'res_name':res_name_test,
                        'atom_ids':res_atom_ids,
                        'atom_names':atom_names_list,
                        'elements':elements_list,
                        'coords':coords_list
                        })            


    n_atoms=len(res_atom_ids)
    # finds bonds involved only with atoms inside each residue
    for b in bond_data:
        if b['a_i'] in res_atom_ids and b['a_j'] in res_atom_ids:
            res_bonds.append((b['a_j'], b['a_i']))

    # if there are no atoms found in a residue
    if n_atoms==0:
        warn=f'res id {n} has no atoms'
        warnings.append(warn)
        print(warn)

    # create graph of residue    
    if n_atoms != 0:
        matching=False
        res_name=res_dict[0]['res_name']
        for res in ref_dict:
            ref_name=res['res_name']
            if res_name.startswith(f'{ref_name}'):
                matched_res=ref_name
                graph_res, dict_res = create_graph(res_atom_ids, elements_list, res_bonds, res_weights)
                graph_ref, dict_ref=create_graph(res['nodes'], res['names'], res['edges'], res['weights'])
                matching=True
                break  # Stop searching
        
        if matching:
            matches_raw=[]
            ISMAGS = nx.isomorphism.ISMAGS(graph_res, graph_ref) # importa qué va primero
            LCIS = list(ISMAGS.largest_common_subgraph(symmetry=True))
            
            # Filter to also match chemical elements
            for d in LCIS:
                ele_ref=[]
                ele_res=[]
                for key, value in d.items():
                    ele_ref.append(next((item_ref['name'] for item_ref in dict_ref if item_ref['id'] == value), None))
                    ele_res.append(next((item['element'] for item in res_dict if item['id'] == key), None))
                if ele_ref == ele_res:
                    matches_raw.append(d)
            
            # searching non-overlapping solutions
            u_matches=find_uniques_matches(matches_raw, len(graph_res.nodes()))

            if len(u_matches) > 1:
                warn=f'Warning: {len(u_matches)} matches found for {res_name} {n}'
                warnings.append(warn)
                print(warn)

            if len(u_matches) == 0:
                warn=f'Warning: no matches found for {res_name} {n}'
                warnings.append(warn)
                print(warn)
            
            
            # appends unique solutions found (this relates atom id with reference dictionary)
            for u in u_matches:
                matches_in_total.append({'res_name':matched_res, 'dictionary':u, 'res_id':n})

        
        if not matching:
            warn=f'residue name {res_name} with res_id {n} was not found in reference dictionary'
            warnings.append(warn)
            print(warn)


# find bonds only within the pocket
pocket_bonds=[]
for b in bond_data:
    if b['a_i'] in pocket_ids and b['a_j'] in pocket_ids:        
        pocket_bonds.append({
            'a_i': b['a_i'],
            'a_j': b['a_j']})


# ________________ original to reference dictionary of pocket ________________ 
final_dict=[]
for m in matches_in_total:
    for ref in ref_dict:
        
        # for each residue
        if m['res_name']==ref['res_name']:            
            atom_id_list=[]
            type_list=[]
            coord_list=[]

            # translation indices (v is ref atom id, k is original id)
            for v, k in zip(m['dictionary'][0].values(),
                            m['dictionary'][0].keys()):

                for node,name in zip(ref['nodes'],ref['type']):
                    if node == v:
                        atom_id_list.append(k)
                        type_list.append(name)

                for a in atom_data:
                    if k == a['id']:
                        coord_list.append(a['coords'])
            
            final_dict.append({'res_id':m['res_id'],
                               'res_name':ref['res_name'],
                               'atom_ids':atom_id_list,
                               'type':type_list,
                               'coords':coord_list, 
                               })

# ___________________________ Finding non bonded residues ___________________________

N_id_list=[] # id of N atoms inside pocket connected outside of pocket
C_id_list=[] # id of C atoms inside pocket connected outside of pocket

# selection of non bonded residues
N_nonbond_res=[] # list of residues not connected by N terminus
C_nonbond_res=[] # list of residues not connected by C terminus

N_bonding=[] # ids of outise pocket atoms N Connected
C_bonding=[] # ids of outise pocket atoms C Connected

for d in final_dict:
    atoms=d['atom_ids']
    for i in range(len(atoms)):

        # searching non connected N residues
        if d['type'][i]=='N':
            N_id=d['atom_ids'][i]
            for b in bond_data:
                # search only for bonds related with name N
                if (b['a_i'] == N_id or b['a_j'] == N_id):
                    # search only for bonds wit atoms outside of pocket
                    if b['a_i'] not in pocket_ids or b['a_j'] not in pocket_ids:
                        N_nonbond_res.append(d['res_id'])
                        N_bonding.append(b['a_i'])
                        N_bonding.append(b['a_j'])
                        N_bonding.remove(N_id)
                        N_id_list.append(N_id)
                        
        # searching non connected C residues
        if d['type'][i]=='C':
            C_id=d['atom_ids'][i]
            for b in bond_data:
                # search only for bonds related with name N
                if (b['a_i'] == C_id or b['a_j'] == C_id):
                    # search only for bonds wit atoms outside of pocket
                    if b['a_i'] not in pocket_ids or b['a_j'] not in pocket_ids:
                        C_nonbond_res.append(d['res_id'])
                        C_bonding.append(b['a_i'])
                        C_bonding.append(b['a_j'])
                        C_bonding.remove(C_id)
                        C_id_list.append(C_id)


# __________________ Coordinates search _____________________
# carbon bonding
CH_coords=[]
for C_p_id, C_b_id in zip(C_id_list, C_bonding):
    for a in atom_data:
        if a['id']==C_p_id:
            C_p_coord=a['coords']
        if a['id']==C_b_id:
            C_b_coord=a['coords']
    bond_vector=C_b_coord-C_p_coord
    unit_vector = bond_vector/np.linalg.norm(bond_vector)
    H_coord=unit_vector*r_CH+C_p_coord
    CH_coords.append(np.round(H_coord, decimals=3))

# Nitrogen bonding
NH_coords=[]
for N_p_id, N_b_id in zip(N_id_list, N_bonding):
    for a in atom_data:
        if a['id']==N_p_id:
            N_p_coord=a['coords']
        if a['id']==N_b_id:
            N_b_coord=a['coords']
    bond_vector=N_b_coord-N_p_coord
    unit_vector = bond_vector/np.linalg.norm(bond_vector)
    H_coord=unit_vector*r_NH+N_p_coord
    NH_coords.append(np.round(H_coord, decimals=3))

# __________________ write coordinates __________________

# mol2 format
c_01=1 # id of fisrt atom
c_02=0 # counter for HC atoms
c_03=0 # counter for HN atoms

lines_to_write=[]
# ligand atoms
for p in lig_dict:
    r_l=1 # ligand is first residue
    for i in range(len(p['atom_ids'])):
        line=mol_atom_writer(c_01, p['atom_names'][i],p['coords'][i],
                        p['elements'][i], r_l, p['res_name'])
        lines_to_write.append(line)
        c_01+=1

# pocket atoms
for r_i, p in enumerate(pocket_dict):
    r_n=r_l+r_i+1 # residue id
    for i in range(len(p['atom_ids'])):
        line=mol_atom_writer(c_01, p['atom_names'][i],p['coords'][i],
                        p['elements'][i], r_n, p['res_name'])
        lines_to_write.append(line)
        c_01+=1
    # Carbon hydrogens
    if p['res_id'] in C_nonbond_res:
        line=mol_atom_writer(c_01, 'HC', CH_coords[c_02],
                        'H', r_n, p['res_name'])
        lines_to_write.append(line)
        c_01+=1
        c_02+=1

    if p['res_id'] in N_nonbond_res:
        line=mol_atom_writer(c_01, 'HN', NH_coords[c_03],
                        'H', r_n, p['res_name'])
        lines_to_write.append(line)

        c_01+=1
        c_03+=1


last_res_id=r_n
n_atoms=len(lines_to_write)

# write mol2 format
write_mol2(output_file, mol_name, n_atoms, last_res_id, lines_to_write)

# write warnings file
with open('warnings.txt', 'w') as f:
        for w in warnings:
            f.write(w + '\n')
