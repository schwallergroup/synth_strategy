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
    Detects if the synthesis route uses a carbamate-based protection strategy (e.g., Cbz, Boc, Fmoc).
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

    # Track protection/deprotection events
    protection_events = []
    deprotection_events = []
    cbz_protected_molecules = []
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if molecule has Cbz group (carbamic ester)
            if checker.check_fg("Carbamic ester", mol_smiles):
                cbz_protected_molecules.append((mol_smiles, depth))
                if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check for Cbz protection (adding Cbz to amine)
            has_cbz_in_reactants = any(checker.check_fg("Carbamic ester", r) for r in reactants)
            has_cbz_in_product = checker.check_fg("Carbamic ester", product_part)

            if not has_cbz_in_reactants and has_cbz_in_product:
                protection_events.append((rsmi, depth))
                if "protection_Carbamic ester" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("protection_Carbamic ester")

            # Check for Cbz deprotection (removing Cbz from amine)
            if has_cbz_in_reactants and not has_cbz_in_product:
                deprotection_events.append((rsmi, depth))
                if "deprotection_Carbamic ester" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("deprotection_Carbamic ester")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Criteria for Cbz protection strategy:
    # 1. At least one protection event at high depth (early in synthesis)
    # 2. At least one deprotection event at low depth (late in synthesis)
    # 3. Cbz group maintained through multiple steps (at least 3 depths)

    if protection_events and deprotection_events:
        # Find earliest protection (highest depth) and latest deprotection (lowest depth)
        earliest_protection_depth = max([depth for _, depth in protection_events])
        latest_deprotection_depth = min([depth for _, depth in deprotection_events])

        # Check if protection happens early and deprotection happens late
        if (
            earliest_protection_depth > max_depth * 0.6
            and latest_deprotection_depth < max_depth * 0.4
        ):
            # Check if Cbz is maintained through multiple steps
            cbz_depths = sorted([depth for _, depth in cbz_protected_molecules])
            if len(cbz_depths) >= 3:
                result = True
                # Add structural constraints if all conditions met
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["protection_Carbamic ester", "deprotection_Carbamic ester"]}})
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "protection_Carbamic ester", "position": "early_stage"}})
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "deprotection_Carbamic ester", "position": "late_stage"}})
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Carbamic ester", "operator": ">=", "value": 3}})

    # Alternative check: if we have Cbz groups maintained through multiple steps (at least 3)
    # even without clear protection/deprotection events
    if not result and len(cbz_protected_molecules) >= 3:
        cbz_depths = sorted([depth for _, depth in cbz_protected_molecules])
        depth_span = max(cbz_depths) - min(cbz_depths)
        if depth_span >= 2:  # Maintained through at least 3 consecutive depths
            result = True
            # Add structural constraint for count if this path triggers the result
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Carbamic ester", "operator": ">=", "value": 3}})

    return result, findings_json