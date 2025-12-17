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
    This function detects late-stage alcohol mesylation (conversion of alcohol to mesylate
    in the final steps of the synthesis)
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

    # Track if mesylation was found and at what depth
    mesylation_found = False
    mesylation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal mesylation_found, mesylation_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Robustness check for invalid SMILES
            if product_mol is None or any(m is None for m in reactant_mols):
                return

            # SMARTS string for alcohol (primary/secondary/phenol)
            alcohol_pattern = "[#6][#8;H1]"

            # SMARTS string for mesylate group
            mesylate_pattern = "[#6][#8][#16](=[#8])(=[#8])[#6]"

            # Check for the specific conversion of an alcohol to a mesylate
            is_alcohol_consumed = checker.is_group_consumed(product_mol, reactant_mols, alcohol_pattern)
            is_mesylate_formed = checker.is_group_formed(product_mol, reactant_mols, mesylate_pattern)

            if is_alcohol_consumed:
                if "alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("alcohol")
            if is_mesylate_formed:
                if "mesylate" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("mesylate")

            if is_alcohol_consumed and is_mesylate_formed:
                mesylation_found = True
                mesylation_depth = min(mesylation_depth, depth)
                if "alcohol_mesylation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("alcohol_mesylation")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if mesylation was found in the late stage (depth <= 1)
    late_stage = mesylation_found and mesylation_depth <= 1

    if late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "alcohol_mesylation",
                "position": "late_stage"
            }
        })

    print(f"Late-stage alcohol mesylation detected: {late_stage}")
    print(f"Mesylation depth: {mesylation_depth if mesylation_found else 'Not found'}")

    return late_stage, findings_json
