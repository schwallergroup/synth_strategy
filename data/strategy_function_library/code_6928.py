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
    Detects if the synthetic route involves transformation of a carboxylic acid
    to a ketone during ring formation.
    """
    found_transformation = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal found_transformation, findings_json

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check product for ketone
                if checker.check_fg("Ketone", product_smiles):
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    product_ring_info = product_mol.GetRingInfo()
                    product_ring_count = product_ring_info.NumRings()

                    # Check each reactant
                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg(
                            "Carboxylic acid", reactant_smiles
                        ):
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            # Check for negation constraint: Ketone not in reactant with carboxylic acid
                            if not checker.check_fg("Ketone", reactant_smiles):
                                # This condition satisfies the 'negation' structural constraint
                                # { "type": "negation", "details": { "target": "Ketone", "scope": "reactant_with_carboxylic_acid" }}
                                # We don't add this to findings_json directly as it's a condition for the main finding.
                                pass
                            else:
                                # If Ketone IS in the reactant, this specific path doesn't satisfy the negation constraint
                                continue # Skip this reactant as it doesn't fit the desired pattern

                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            reactant_ring_info = reactant_mol.GetRingInfo()
                            reactant_ring_count = reactant_ring_info.NumRings()

                            # Check if product has more rings than reactant (ring formation)
                            if product_ring_count > reactant_ring_count:
                                # At this point, we have confirmed a carboxylic acid-containing
                                # reactant forms a ketone-containing product with an increase
                                # in ring count. This is the target transformation.
                                found_transformation = True
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "targets": [
                                            "Carboxylic acid",
                                            "Ketone",
                                            "ring_formation"
                                        ]
                                    }
                                })
                                findings_json["structural_constraints"].append({
                                    "type": "negation",
                                    "details": {
                                        "target": "Ketone",
                                        "scope": "reactant_with_carboxylic_acid"
                                    }
                                })

                                print(
                                    f"Found carboxylic acid to ketone transformation with ring formation at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                print(f"Reactant with carboxylic acid: {reactant_smiles}")
                                print(f"Product with ketone: {product_smiles}")
                                return  # Return early once we find a match
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)
            if found_transformation:
                return  # Stop traversal once we find a match

    # Start traversal
    dfs_traverse(route)

    print(f"Carboxylic acid to ketone transformation detected: {found_transformation}")
    return found_transformation, findings_json
