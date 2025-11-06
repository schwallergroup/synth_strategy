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


PRESERVED_FGS_OF_INTEREST = ["Nitrile", "Ester", "Carboxylic acid", "Amide", "Alcohol"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving ring opening in the late stage of synthesis
    while preserving key functional groups.
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

    # Track if ring opening occurs in late stage
    late_stage_ring_opening = False

    # Track if key functional groups are preserved
    key_functional_groups_preserved = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ring_opening, key_functional_groups_preserved, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for ring opening in late stage (depth ≤ 3)
                if depth <= 3 and not late_stage_ring_opening:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        product_ring_info = product_mol.GetRingInfo()
                        product_ring_count = product_ring_info.NumRings()

                        # Check each reactant for more rings than the product
                        for reactant_smiles in reactants_smiles:
                            if reactant_smiles:
                                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                                if reactant_mol:
                                    reactant_ring_info = reactant_mol.GetRingInfo()
                                    reactant_ring_count = reactant_ring_info.NumRings()

                                    if reactant_ring_count > product_ring_count:
                                        # Ring destruction detected
                                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                                        # Positional constraint: late_stage
                                        if depth <= 3:
                                            findings_json["structural_constraints"].append({
                                                "type": "positional",
                                                "details": {
                                                    "target": "ring_destruction",
                                                    "position": "late_stage",
                                                    "condition": "depth <= 3",
                                                    "note": "The ring destruction event must occur within the last 3 steps of the synthesis."
                                                }
                                            })

                                        # Check if any key functional group is preserved
                                        preserved_groups = []
                                        for fg in PRESERVED_FGS_OF_INTEREST:
                                            if checker.check_fg(
                                                fg, product_smiles
                                            ) and checker.check_fg(fg, reactant_smiles):
                                                preserved_groups.append(fg)
                                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                                    findings_json["atomic_checks"]["functional_groups"].append(fg)

                                        if preserved_groups:
                                            key_functional_groups_preserved = True
                                            late_stage_ring_opening = True
                                            # Co-occurrence constraint: ring_destruction and functional_group_preservation
                                            findings_json["structural_constraints"].append({
                                                "type": "co-occurrence",
                                                "details": {
                                                    "targets": [
                                                        "ring_destruction",
                                                        "functional_group_preservation"
                                                    ],
                                                    "scope": "same_reaction",
                                                    "note": "Checks if a ring destruction event occurs in the same reaction where at least one key functional group is preserved (present in both a reactant and the product)."
                                                }
                                            })
                                            break

            except Exception as e:
                # Silently ignore errors in processing a single reaction
                pass

        # Process children
        for child in node.get("children", []):
            # Determine new_depth based on the current node's type
            if node["type"] == "reaction":
                new_depth = depth
            else:  # chemical node
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = late_stage_ring_opening and key_functional_groups_preserved

    return strategy_present, findings_json