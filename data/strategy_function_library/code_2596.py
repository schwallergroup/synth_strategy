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
    Detects if a Boc-protected amine is maintained throughout the synthesis.
    This means the Boc group is present in the final product and is not added
    or removed during the synthesis.
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

    result = True

    # Check if target molecule has Boc group
    if route["type"] == "mol":
        target_mol = route["smiles"]
        target_has_boc = checker.check_fg("Boc", target_mol)
        if not target_has_boc:
            result = False
        else:
            findings_json["atomic_checks"]["functional_groups"].append("Boc")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Boc",
                    "position": "final_product"
                }
            })

    # Track if Boc is maintained in all reactions
    boc_maintained = True

    def dfs_traverse(node, depth=0):
        nonlocal boc_maintained, result, findings_json
        if not boc_maintained:
            return

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_has_boc = checker.check_fg("Boc", product)
            reactant_has_boc = any(checker.check_fg("Boc", reactant) for reactant in reactants)

            # The "maintained" strategy is violated if the Boc group's presence
            # changes between reactants and product.
            if product_has_boc != reactant_has_boc:
                boc_maintained = False
                result = False
                # If Boc is added (reactant_has_boc is False, product_has_boc is True)
                if not reactant_has_boc and product_has_boc:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc protection")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Boc protection"
                        }
                    })
                # If Boc is removed (reactant_has_boc is True, product_has_boc is False)
                elif reactant_has_boc and not product_has_boc:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc deprotection")
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Boc deprotection"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if not boc_maintained:
        result = False

    return result, findings_json