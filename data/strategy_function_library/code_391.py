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
    Detects if a methoxy group is preserved throughout the entire synthesis route.
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

    # Track if we've found a methoxy group that's preserved
    target_has_methoxy = False
    all_methoxy_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal target_has_methoxy, all_methoxy_preserved, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule has a methoxy group
            has_methoxy = False
            if checker.check_fg("Ether", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Look for -OCH3 pattern (methoxy group)
                    methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                    if mol.HasSubstructMatch(methoxy_pattern):
                        has_methoxy = True
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")

            if depth == 0:
                if has_methoxy:
                    target_has_methoxy = True
                    # Structural constraint: positional - target has Ether
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Ether",
                            "position": "last_stage"
                        }
                    })
                else:
                    # If target doesn't have methoxy, no need to check preservation
                    return

        elif node["type"] == "reaction" and target_has_methoxy:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check if product has methoxy
            product_has_methoxy = False
            if checker.check_fg("Ether", product):
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                    if product_mol.HasSubstructMatch(methoxy_pattern):
                        product_has_methoxy = True
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")

            # Check if at least one reactant has methoxy
            reactant_has_methoxy = False
            for r in reactants:
                if checker.check_fg("Ether", r):
                    reactant_mol = Chem.MolFromSmiles(r)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                    ):
                        reactant_has_methoxy = True
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")
                        break

            if product_has_methoxy and not reactant_has_methoxy:
                all_methoxy_preserved = False
                # Structural constraint: negation - de_novo_formation of Ether
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "Ether",
                        "event": "de_novo_formation"
                    }
                })

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
            new_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Only return True if target has methoxy and all methoxy groups are preserved
    result = target_has_methoxy and all_methoxy_preserved
    return result, findings_json
