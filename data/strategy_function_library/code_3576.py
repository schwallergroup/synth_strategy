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
    Detects if the synthesis route uses a nitrile as a key intermediate,
    particularly in the penultimate step.
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

    result = False

    nitrile_nodes = []
    nitrile_reactions = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal result, findings_json
        if path is None:
            path = []

        # For molecule nodes, check if they contain a nitrile group
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            is_nitrile = checker.check_fg("Nitrile", mol_smiles)

            if is_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                # Store node info for nitrile-containing molecules
                nitrile_nodes.append(
                    {
                        "smiles": mol_smiles,
                        "depth": depth,
                        "in_stock": node.get("in_stock", False),
                        "path_index": len(path),
                    }
                )
                print(f"Found nitrile at depth {depth}: {mol_smiles}")

        # For reaction nodes, check if they involve nitrile transformation
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains a nitrile
            reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
            product_has_nitrile = checker.check_fg("Nitrile", product)

            if reactant_has_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
            if product_has_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            # Check for reactions that commonly use nitriles
            nitrile_to_amide = checker.check_reaction("Nitrile to amide", rsmi)
            if nitrile_to_amide:
                findings_json["atomic_checks"]["named_reactions"].append("Nitrile to amide")

            nitrile_to_carboxylic_acid = checker.check_reaction(
                "Oxidation of nitrile to carboxylic acid", rsmi
            )
            if nitrile_to_carboxylic_acid:
                findings_json["atomic_checks"]["named_reactions"].append("Oxidation of nitrile to carboxylic acid")

            if (
                nitrile_to_amide
                or nitrile_to_carboxylic_acid
            ):
                nitrile_reactions.append(
                    {
                        "rsmi": rsmi,
                        "depth": depth,
                        "reactant_has_nitrile": reactant_has_nitrile,
                        "product_has_nitrile": product_has_nitrile,
                    }
                )
                print(f"Found nitrile-related reaction at depth {depth}")

        # Continue traversing the tree
        new_path = path + [node]
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            if node["type"] == "reaction":
                # Depth remains the same when going from reaction to chemical
                next_depth = depth
            else:
                # Depth increases when going from chemical to reaction
                next_depth = depth + 1
            dfs_traverse(child, next_depth, new_path)

    dfs_traverse(route)

    # Check for strategic use of nitrile intermediates
    if not nitrile_nodes:
        print("No nitrile-containing molecules found")
        return False, findings_json

    # Filter out nitriles that are starting materials
    intermediate_nitriles = [n for n in nitrile_nodes if not n.get("in_stock", False)]
    if not intermediate_nitriles:
        print("Found nitriles, but they are all starting materials")
        findings_json["structural_constraints"].append({"type": "negation", "details": {"event": "all_nitriles_are_starting_materials"}})
        return False, findings_json

    # Check if any nitrile is used in a reaction (consumed as reactant)
    for reaction in nitrile_reactions:
        if reaction["reactant_has_nitrile"] and not reaction["product_has_nitrile"]:
            print("Found a reaction where nitrile is consumed")
            findings_json["atomic_checks"]["named_reactions"].append("nitrile_consumption")
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "nitrile_consumption", "operator": ">=", "value": 1}})
            result = True
            return result, findings_json

    # Check specifically for nitrile at penultimate step (depth 1)
    penultimate_nitriles = [n for n in intermediate_nitriles if n["depth"] == 1]
    if penultimate_nitriles:
        print("Found nitrile at penultimate step")
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitrile", "context": "intermediate", "position": "penultimate_stage"}})
        result = True
        return result, findings_json

    # Check for any intermediate nitriles in the middle of the synthesis
    middle_nitriles = [n for n in intermediate_nitriles if 0 < n["depth"] < 3]
    if middle_nitriles:
        print("Found nitrile as intermediate in the middle of synthesis")
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitrile", "context": "intermediate", "position": "middle_stages"}})
        result = True
        return result, findings_json

    print("Nitriles found but not used as key intermediates")
    return result, findings_json
