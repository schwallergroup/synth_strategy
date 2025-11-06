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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving formylation in the late stage of synthesis
    (low depth in the retrosynthetic tree).
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

    formylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal formylation_detected, findings_json

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions (depth 0 or 1)
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Extract reactants and product
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for formyl group in product
                product_has_formyl = checker.check_fg(
                    "Aldehyde", product_smiles
                )
                if product_has_formyl:
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

                if product_has_formyl:
                    print(f"Product has formyl group: {product_smiles}")

                    # Check if main reactant doesn't have formyl group
                    main_reactant_has_formyl = False
                    for reactant in reactants_smiles:
                        if reactant and Chem.MolFromSmiles(reactant):
                            # Skip small formylating agents
                            if len(reactant) > 5 and checker.check_fg("Aldehyde", reactant):
                                main_reactant_has_formyl = True
                                if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                                break

                    # Check for common formylating agents
                    has_formylating_agent = False

                    for r in reactants_smiles:
                        if not r:
                            continue
                        r_mol = Chem.MolFromSmiles(r)
                        if not r_mol:
                            continue

                        # Check if reactant is a small molecule with formyl group
                        if len(r) <= 5 and checker.check_fg("Aldehyde", r):
                            has_formylating_agent = True
                            if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                            print(f"Found formylating agent: {r}")
                            break

                    # Formylation detected if:
                    # 1. Product has formyl group (checked above)
                    # 2. Main reactant doesn't have formyl AND a formylating agent is present
                    if not main_reactant_has_formyl and has_formylating_agent:
                        print(f"Late-stage formylation detected at depth {depth}")
                        formylation_detected = True
                        if "formylation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("formylation")
                        # Add the structural constraint for late_stage formylation
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "formylation",
                                "position": "late_stage"
                            }
                        })

        # Traverse children with increased depth
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return formylation_detected, findings_json
