import rdkit.Chem as Chem
from rdkit.Chem import AllChem
class Check:
    def __init__(self, fg_dict, reaction_dict, ring_dict):
        self.fg_dict = fg_dict
        self.reaction_class_dict = reaction_dict
        self.ring_dict = ring_dict

    def check_ring(self, name, smi):
        ring_smiles = self.ring_dict[name]
        if ring_smiles == []:
            print(f"No smiles found for {name}")
            return False
        for ring_smi in ring_smiles:
            if Chem.MolFromSmiles(smi).HasSubstructMatch(Chem.MolFromSmiles(ring_smi)):
                return True
        return False
    
    def get_ring_smiles(self, name):
        ring_smiles = self.ring_dict[name]
        if ring_smiles == []:
            print(f"No smiles found for {name}")
            return []
        return ring_smiles[0]
    
    def get_ring_atom_indices(self, ring_name, query_mol_smiles):
        """
        Find the atom indices in the query molecule where a specific ring pattern matches.
        
        Parameters:
        ring_name (str): The name of the ring to search for
        query_mol_smiles (str): SMILES string of the molecule to search in
        
        Returns:
        list: List of lists containing atom indices for each ring match, or empty list if no matches
        """
        ring_smiles_list = self.ring_dict[ring_name]
        if not ring_smiles_list:
            print(f"No SMILES patterns found for {ring_name}")
            return []
        
        # Create molecule from SMILES
        query_mol = Chem.MolFromSmiles(query_mol_smiles)
        if query_mol is None:
            print(f"Could not create molecule from SMILES: {query_mol_smiles}")
            return []
        
        all_matches = []
        for ring_smi in ring_smiles_list:
            # Create query molecule from ring SMILES
            ring_query = Chem.MolFromSmiles(ring_smi)
            if ring_query is None:
                print(f"Could not create query from SMILES: {ring_smi}")
                continue
            
            # Find all matches
            matches = query_mol.GetSubstructMatches(ring_query)
            if matches:
                all_matches.extend(list(matches))
        
        return all_matches


    def check_fg(self, fg, smi):
        fg_smarts_list = self.fg_dict[fg]
        if fg_smarts_list == []:
            print(f"No smarts found for {fg}")
            return False
        for fg_smarts in fg_smarts_list:
            if Chem.MolFromSmiles(smi).HasSubstructMatch(Chem.MolFromSmarts(fg_smarts)):
                return True
        return False

    def check_reaction(self, rxn_type, rsmi):
        rxn_smarts_list = self.reaction_class_dict[rxn_type]
        if rxn_smarts_list == []:
            print(f"No smarts found for {rxn_type}")
            return False
        for rxn_smarts in rxn_smarts_list:
            rxn_smi = AllChem.ReactionFromSmarts(rsmi, {}, True)
            rxn_smart_obj = AllChem.ReactionFromSmarts(rxn_smarts)
            if AllChem.HasReactionSubstructMatch(rxn_smart_obj, rxn_smi):
                return True
        return False
    
    def get_fg_smart(self, fg):
        return self.fg_dict[fg]
    
    def get_reaction_smart(self, rxn_type):
        return self.reaction_class_dict[rxn_type]
    
    def get_fg_atom_indices(self, fg_name, query_mol_smiles):
        """
        Find the atom indices in the query molecule where a specific functional group matches.
        
        Parameters:
        fg_name (str): The name of the functional group to search for
        query_mol_smiles (str): SMILES string of the molecule to search in
        
        Returns:
        list: List of lists containing atom indices for each match, or empty list if no matches
        """
        fg_smarts_list = self.fg_dict[fg_name]
        if not fg_smarts_list:
            print(f"No SMARTS patterns found for {fg_name}")
            return []
        
        # Create molecule from SMILES
        query_mol = Chem.MolFromSmiles(query_mol_smiles)
        if query_mol is None:
            print(f"Could not create molecule from SMILES: {query_mol_smiles}")
            return []
        
        all_matches = []
        
        for fg_smarts in fg_smarts_list:
            fg_query = Chem.MolFromSmarts(fg_smarts)
            if fg_query is None:
                print(f"Could not create query from SMARTS: {fg_smarts}")
                continue
            
            # Find all matches
            matches = query_mol.GetSubstructMatches(fg_query)
            if matches:
                all_matches.extend(list(matches))
        
        return all_matches