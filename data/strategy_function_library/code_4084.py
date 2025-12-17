from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of a new bond between an aromatic ring and an sp³-hybridized carbon atom.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    aryl_sp3_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal aryl_sp3_formation_found, findings_json

        if aryl_sp3_formation_found:
            return  # Early return if already found

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            try:
                reactants_part, product_part = rsmi.split(">>", 1)

                if not product_part or not reactants_part:
                    return

                product_mol = Chem.MolFromSmiles(product_part)
                if not product_mol:
                    return

                # Find all C(sp³)-Aryl connections in the product by atom map numbers
                aryl_sp3_connections = []
                for atom in product_mol.GetAtoms():
                    if (
                        atom.GetSymbol() == "C"
                        and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3
                    ):
                        atom_map_num = atom.GetAtomMapNum()
                        if atom_map_num > 0:
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetIsAromatic():
                                    neighbor_map_num = neighbor.GetAtomMapNum()
                                    if neighbor_map_num > 0:
                                        aryl_sp3_connections.append(
                                            (atom_map_num, neighbor_map_num)
                                        )

                if not aryl_sp3_connections:
                    return

                # Check if these connections exist in any of the reactants
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                reactant_mols = [m for m in reactant_mols if m]

                for sp3_map, aryl_map in aryl_sp3_connections:
                    connection_exists_in_reactants = False
                    for reactant_mol in reactant_mols:
                        sp3_atom = None
                        aryl_atom = None
                        for atom in reactant_mol.GetAtoms():
                            if atom.GetAtomMapNum() == sp3_map:
                                sp3_atom = atom
                            elif atom.GetAtomMapNum() == aryl_map:
                                aryl_atom = atom
                        
                        if sp3_atom and aryl_atom:
                            # Check if they are bonded in this reactant
                            if reactant_mol.GetBondBetweenAtoms(sp3_atom.GetIdx(), aryl_atom.GetIdx()) is not None:
                                connection_exists_in_reactants = True
                                break # Found the bond, no need to check other reactants
                    
                    if not connection_exists_in_reactants:
                        aryl_sp3_formation_found = True
                        findings_json["atomic_checks"]["named_reactions"].append("aryl_sp3_bond_formation")
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "aryl_sp3_bond_formation",
                                "operator": ">=",
                                "value": 1
                            }
                        })
                        return # Found a newly formed bond, can exit

            except Exception:
                # Error during RDKit processing, skip this reaction
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return aryl_sp3_formation_found, findings_json
