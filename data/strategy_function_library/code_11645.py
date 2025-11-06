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


UREA_LIKE_FGS = ["Urea", "Thiourea"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects late-stage formation of specific urea-like functional groups, including Urea and Thiourea, that connects two complex fragments.
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

    urea_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected, findings_json

        if node["type"] == "reaction" and depth <= 3:  # Late stage
            # Record late_stage constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "functional_group_formation",
                    "position": "late_stage",
                    "condition": "depth <= 3"
                }
            })
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains a target FG
                    has_urea_product = False
                    for fg in UREA_LIKE_FGS:
                        if checker.check_fg(fg, product):
                            has_urea_product = True
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                    if has_urea_product:
                        # Check if target FG is formed (not present in reactants)
                        urea_in_reactants = False
                        for r in reactants:
                            for fg in UREA_LIKE_FGS:
                                if checker.check_fg(fg, r):
                                    urea_in_reactants = True
                                    break
                            if urea_in_reactants:
                                break

                        if not urea_in_reactants:
                            # Record functional_group_formation if it's a new formation
                            findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")

                            # Check if this reaction connects two complex fragments
                            complex_reactants = 0
                            for reactant in reactants:
                                try:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        # Heuristic for complexity: >5 atoms or has rings
                                        if (
                                            reactant_mol.GetNumAtoms() > 5
                                            or reactant_mol.GetRingInfo().NumRings() > 0
                                        ):
                                            complex_reactants += 1
                                except Exception:
                                    pass  # Ignore SMILES parsing errors

                            # The strategy requires connecting at least two complex fragments.
                            if complex_reactants >= 2:
                                # Record complex_reactants constraint if met
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "complex_reactants",
                                        "operator": ">=",
                                        "value": 2
                                    }
                                })
                                urea_formation_detected = True
            except Exception:
                # Broad exception to handle any unexpected errors during analysis
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return urea_formation_detected, findings_json