from typing import Tuple, Dict, List
import copy
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
    Detects if the synthesis maintains an aromatic ether throughout the route.

    This strategy checks if:
    1. The final product contains an aromatic ether
    2. The aromatic ether is maintained through all reaction steps
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

    # Track if we've found the final product with aromatic ether
    final_product_has_aromatic_ether = False

    # Track if aromatic ether is maintained through reactions
    ether_maintained = True

    # Keep track of molecules we've checked
    molecules_checked = 0

    def has_aromatic_ether(smiles):
        """Check if a molecule has an aromatic ether group"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Check if molecule has ether functional group
            if checker.check_fg("Ether", smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Ether")
                # Check if the ether is connected to an aromatic carbon
                aromatic_carbon_pattern = Chem.MolFromSmarts("c-[O]")
                return mol.HasSubstructMatch(aromatic_carbon_pattern)
        return False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_aromatic_ether, ether_maintained, molecules_checked, findings_json

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            molecules_checked += 1

            # Check if this is the final product (depth 0)
            if depth == 0:
                final_product_has_aromatic_ether = has_aromatic_ether(mol_smiles)
                if final_product_has_aromatic_ether:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "aromatic_ether", "position": "last_stage"}})
                print(
                    f"Final product {'has' if final_product_has_aromatic_ether else 'does not have'} aromatic ether: {mol_smiles}"
                )

            # For reaction nodes, check if aromatic ether is maintained
            for child in node.get("children", []):
                if child["type"] == "reaction" and child.get("metadata", {}).get("rsmi"):
                    rsmi = child["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if product has aromatic ether
                    product_has_ether = has_aromatic_ether(product)

                    # Check if any reactant has aromatic ether
                    reactant_has_ether = any(has_aromatic_ether(r) for r in reactants)

                    # If a reactant has an ether but the product does not, it was cleaved.
                    if reactant_has_ether and not product_has_ether:
                        print(f"Aromatic ether not maintained in reaction: {rsmi}")
                        ether_maintained = False
                        # Assuming 'aromatic_ether_cleavage' is a named reaction that would be detected here
                        # This is a conceptual link, as the original code doesn't explicitly check for a reaction name.
                        # We'll add it to findings if the condition for cleavage is met.
                        if "aromatic_ether_cleavage" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatic_ether_cleavage")
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "aromatic_ether_cleavage"}})

        # Continue traversal
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True only if:
    # 1. We checked at least one molecule
    # 2. The final product has aromatic ether
    # 3. The aromatic ether is maintained through all reactions
    result = False
    if molecules_checked > 0 and final_product_has_aromatic_ether and ether_maintained:
        print("Aromatic ether maintained throughout synthesis")
        result = True

    return result, findings_json
