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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis involves a Boc protection/deprotection sequence for an amine by identifying specific named reactions and verifying they occur in a logical sequence on the same molecular scaffold. It checks for protection reactions from the list: ['Boc amine protection', 'Boc amine protection explicit', 'Boc amine protection with Boc anhydride', 'Boc amine protection (ethyl Boc)', 'Boc amine protection of secondary amine', 'Boc amine protection of primary amine']. It checks for deprotection reactions from the list: ['Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', 'Tert-butyl deprotection of amine'].
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

    # Track protection and deprotection events with their depths
    protection_events = []  # Will store (depth, atom_maps) tuples
    deprotection_events = []  # Will store (depth, atom_maps) tuples

    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal result # Declare result as nonlocal to modify it
        nonlocal findings_json # Declare findings_json as nonlocal to modify it

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection reactions
                is_boc_protection = False
                for name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_boc_protection = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                # Check for Boc deprotection reactions
                is_boc_deprotection = False
                for name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_boc_deprotection = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                # Extract atom maps for protection events
                if is_boc_protection:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        atom_maps = []
                        for atom in product_mol.GetAtoms():
                            map_num = atom.GetAtomMapNum()
                            if map_num > 0:
                                atom_maps.append(map_num)
                        if atom_maps:
                            protection_events.append((depth, set(atom_maps)))

                # Extract atom maps for deprotection events
                if is_boc_deprotection:
                    boc_reactant_mol = None
                    for reactant in reactants:
                        if checker.check_fg("Boc", reactant):
                            if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Boc")
                            boc_reactant_mol = Chem.MolFromSmiles(reactant)
                            break

                    if boc_reactant_mol:
                        atom_maps = []
                        for atom in boc_reactant_mol.GetAtoms():
                            map_num = atom.GetAtomMapNum()
                            if map_num > 0:
                                atom_maps.append(map_num)
                        if atom_maps:
                            deprotection_events.append((depth, set(atom_maps)))

            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    if protection_events and deprotection_events:
        # Add co-occurrence constraint if both types of events are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Boc_protection_reaction",
                    "Boc_deprotection_reaction"
                ],
                "description": "The route must contain at least one reaction from the set of defined Boc protection reactions and at least one from the set of defined Boc deprotection reactions."
            }
        })

        # Check if there's a valid protection/deprotection sequence
        # In retrosynthetic direction, deprotection should be at lower depth than protection
        for deprotect_depth, deprotect_maps in deprotection_events:
            for protect_depth, protect_maps in protection_events:
                # Check if depths are in correct order (deprotection before protection in retrosynthesis)
                if deprotect_depth < protect_depth:
                    # Check if there's overlap in atom maps (same molecular fragment)
                    if deprotect_maps.intersection(protect_maps):
                        result = True # Set the main result to True
                        # Add sequence constraint if found
                        findings_json["structural_constraints"].append({
                            "type": "sequence",
                            "details": {
                                "before": "Boc_deprotection_reaction",
                                "after": "Boc_protection_reaction",
                                "description": "A Boc deprotection reaction must occur at a lower retrosynthetic depth (closer to the final product) than a corresponding Boc protection reaction on the same molecular scaffold, verified by overlapping atom map numbers."
                            }
                        })
                        # No need to break, as we want to find all such sequences if possible, but the result is already True

    return result, findings_json
