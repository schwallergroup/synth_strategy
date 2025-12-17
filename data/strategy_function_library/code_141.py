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


BORONIC_ACID_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Chan-Lam alcohol",
    "Chan-Lam amine",
    "Petasis reaction with amines and boronic acids",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of a boronic acid as a key intermediate. This strategy is confirmed if a boronic acid is formed (either from scratch or by deprotecting a boronic ester) and is subsequently consumed in one of a specific list of coupling reactions, including Suzuki, Chan-Lam, and Petasis reactions.
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

    boronic_acid_formed = False
    boronic_acid_used = False

    def dfs_traverse(node, depth=0):
        nonlocal boronic_acid_formed, boronic_acid_used, findings_json

        # Process reaction nodes
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid formation reactions
                product_mol = product
                if checker.check_fg("Boronic acid", product_mol):
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                    print(f"Depth {depth}: Found boronic acid in product: {product_mol}")

                    has_boronic_acid_in_reactants = any(
                        checker.check_fg("Boronic acid", r) for r in reactants
                    )
                    has_boronic_ester_in_reactants = any(
                        checker.check_fg("Boronic ester", r) for r in reactants
                    )

                    if not has_boronic_acid_in_reactants and not has_boronic_ester_in_reactants:
                        boronic_acid_formed = True
                        print(f"Depth {depth}: Confirmed boronic acid formation")
                    elif has_boronic_ester_in_reactants:
                        boronic_acid_formed = True
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                        print(f"Depth {depth}: Confirmed boronic acid formation from ester")

                # Check for boronic acid usage in coupling reactions
                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                        print(f"Depth {depth}: Found boronic acid in reactant: {reactant}")

                        # Check if this is a known coupling reaction
                        for r_name in BORONIC_ACID_COUPLING_REACTIONS:
                            if checker.check_reaction(r_name, rsmi):
                                boronic_acid_used = True
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                print(
                                    f"Depth {depth}: Confirmed boronic acid used in coupling reaction"
                                )
                                break

        # Process children nodes
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    print(
        f"Final result - Boronic acid formed: {boronic_acid_formed}, Boronic acid used: {boronic_acid_used}"
    )
    
    result = boronic_acid_formed and boronic_acid_used
    if result:
        # Add the structural constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "boronic_acid_formation",
                    "boronic_acid_consumption_in_coupling"
                ]
            }
        })

    return result, findings_json