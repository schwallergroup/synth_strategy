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


REDUCTIVE_AMINATION_REACTIONS = [
    "Reductive amination with ketone",
    "Reductive amination with aldehyde",
    "Mignonac reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving the early-stage formation of a chiral amine via specific reductive amination reactions. This strategy is identified if a reaction matches one of the defined `REDUCTIVE_AMINATION_REACTIONS` (e.g., from a ketone or aldehyde, or a Mignonac reaction) and results in the formation of an amine at a chiral center. The reaction must occur in the first half of the synthesis (i.e., `depth > max_depth / 2`).
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

    # Track if we found ketone reduction to amine
    ketone_reduction_found = False
    # Track the depth at which reduction occurs
    reduction_depth = None
    # Track total depth of synthesis
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal ketone_reduction_found, reduction_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node.get("type") == "reaction":
            # Extract reaction information
            if "metadata" in node and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                is_reductive_amination = False
                for reaction_name in REDUCTIVE_AMINATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_reductive_amination = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                amine_fgs = []
                if checker.check_fg("Primary amine", product):
                    amine_fgs.append("Primary amine")
                if checker.check_fg("Secondary amine", product):
                    amine_fgs.append("Secondary amine")
                if checker.check_fg("Tertiary amine", product):
                    amine_fgs.append("Tertiary amine")

                amine_in_product = len(amine_fgs) > 0
                for fg in amine_fgs:
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)

                has_chiral_amine = False
                if is_reductive_amination:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            Chem.AssignStereochemistry(product_mol)
                            chiral_centers = Chem.FindMolChiralCenters(
                                product_mol, includeUnassigned=False
                            )
                            has_chiral_center = len(chiral_centers) > 0

                            if has_chiral_center and amine_in_product:
                                for atom in product_mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() == "N" and atom.GetTotalNumHs() > 0
                                    ):  # Amine nitrogen
                                        for neighbor in atom.GetNeighbors():
                                            neighbor_idx = neighbor.GetIdx()
                                            if any(
                                                center[0] == neighbor_idx for center in chiral_centers
                                            ):
                                                has_chiral_amine = True
                                                break
                                    if has_chiral_amine:
                                        break
                    except Exception:
                        has_chiral_amine = False

                if is_reductive_amination and has_chiral_amine:
                    ketone_reduction_found = True
                    reduction_depth = depth
                    # Add co-occurrence constraint if not already added
                    co_occurrence_constraint = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "reductive_amination",
                                "amine_formation"
                            ]
                        }
                    }
                    if co_occurrence_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(co_occurrence_constraint)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    is_early_stage = reduction_depth is not None and reduction_depth > (max_depth / 2)

    result = False
    if ketone_reduction_found and is_early_stage:
        result = True
        # Add positional constraint if not already added
        positional_constraint = {
            "type": "positional",
            "details": {
                "target": "reductive_amination",
                "position": "first_half"
            }
        }
        if positional_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(positional_constraint)

    return result, findings_json
