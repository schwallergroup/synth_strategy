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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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


COMPLEX_FRAGMENT_FGS = [
    "Carboxylic acid",
    "Ester",
    "Amide",
    "Amine",
    "Alcohol",
    "Ketone",
    "Aldehyde",
    "Nitrile",
    "Halide",
    "Aromatic halide",
]

def main(route, max_depth) -> Tuple[bool, Dict]:
    """
    Identifies key convergent steps by checking for reactions late in the synthesis (depth < 3) that join at least two complex fragments. A fragment is defined as 'complex' if it has more than 10 atoms or more than one ring, and also contains a functional group from the `COMPLEX_FRAGMENT_FGS` list.
    """
    convergent_synthesis_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def has_functional_groups(mol):
        """Check if molecule has significant functional groups"""
        if mol is None:
            return False

        mol_smiles = Chem.MolToSmiles(mol)

        found_fgs = []
        for fg in COMPLEX_FRAGMENT_FGS:
            if checker.check_fg(fg, mol_smiles):
                found_fgs.append(fg)
        
        if found_fgs:
            # Add found functional groups to findings_json, avoiding duplicates
            for fg_name in found_fgs:
                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
            return True

        return False

    def dfs_traverse(reaction, depth, max_depth):
        nonlocal convergent_synthesis_detected, findings_json

        if reaction["type"] == "reaction":
            # Extract reactants for the reaction
            try:
                rsmi = reaction["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if we have at least 2 reactants
                if len(reactants_smiles) >= 2:
                    # Convert to RDKit molecules with error handling
                    reactant_mols = []
                    for r in reactants_smiles:
                        mol = Chem.MolFromSmiles(r)
                        if mol is not None:
                            reactant_mols.append(mol)

                    # Check complexity of each fragment
                    complex_fragments = 0
                    for mol in reactant_mols:
                        # Consider both atom count and ring structures
                        if mol.GetNumAtoms() > 10 or len(Chem.GetSSSR(mol)) > 1:
                            # Check for functional groups
                            if has_functional_groups(mol):
                                complex_fragments += 1

                    # In convergent synthesis, we expect at least 2 complex fragments
                    # at a late stage (low depth)
                    if depth < 3 and complex_fragments >= 2:
                        convergent_synthesis_detected = True
                        reaction_id = reaction.get("metadata", {}).get("ID", "Unknown")
                        print(
                            f"Detected convergent synthesis with {complex_fragments} complex fragments at reaction {reaction_id}, depth {depth}"
                        )
                        # Record structural constraints
                        # Positional constraint: late_stage (depth < 3)
                        if {"type": "positional", "details": {"target": "convergent_reaction", "position": "late_stage (depth < 3)"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "convergent_reaction", "position": "late_stage (depth < 3)"}})
                        # Count constraint: complex_fragments_in_reaction >= 2
                        if {"type": "count", "details": {"target": "complex_fragments_in_reaction", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_fragments_in_reaction", "operator": ">=", "value": 2}})

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in reaction.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if reaction["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same (new_depth = depth)
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal (assuming max_depth is calculated prior to this call)
    # max_depth = calculate_max_depth(route) # Example of what's needed
    dfs_traverse(route, 1, max_depth)

    return convergent_synthesis_detected, findings_json
