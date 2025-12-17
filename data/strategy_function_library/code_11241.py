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


NITROGEN_HETEROCYCLES = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "aziridine", "azetidine", "azepane", "diazepane", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "benzimidazole",
    "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a Boc protection/deprotection strategy on a nitrogen-containing heterocycle. The strategy is identified if a deprotection step is found, along with either a preceding protection step or the presence of a Boc-protected N-heterocycle intermediate. The specific heterocycles checked are defined in the NITROGEN_HETEROCYCLES list.
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

    has_boc_protected_n_heterocycle = False
    has_protection_step = False
    has_deprotection_step = False

    # Track the depth of each node to ensure correct sequence
    protection_depth = float("inf")
    deprotection_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protected_n_heterocycle, has_protection_step, has_deprotection_step
        nonlocal protection_depth, deprotection_depth, findings_json

        if node["type"] == "mol":
            # Check if molecule contains Boc group on a nitrogen heterocycle
            if node.get("smiles"):
                mol_smiles = node["smiles"]

                # Check if molecule has a Boc group
                has_boc = checker.check_fg("Boc", mol_smiles)
                if has_boc:
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")

                # Check if molecule contains a nitrogen heterocycle
                has_n_heterocycle = False
                for ring in NITROGEN_HETEROCYCLES:
                    if checker.check_ring(ring, mol_smiles):
                        has_n_heterocycle = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if has_boc and has_n_heterocycle:
                    has_boc_protected_n_heterocycle = True
                    print(f"Found Boc-protected nitrogen heterocycle: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc protection reactions by name
            protection_reactions = [
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of secondary amine",
                "Boc amine protection of primary amine",
            ]
            for r_name in protection_reactions:
                if checker.check_reaction(r_name, rsmi):
                    has_protection_step = True
                    protection_depth = min(protection_depth, depth)
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Found Boc protection step at depth {depth}: {rsmi}")
                    break

            # Check for Boc deprotection reactions by name
            deprotection_reactions = [
                "Boc amine deprotection",
                "Boc amine deprotection of guanidine",
                "Boc amine deprotection to NH-NH2",
                "Tert-butyl deprotection of amine",
            ]
            for r_name in deprotection_reactions:
                if checker.check_reaction(r_name, rsmi):
                    has_deprotection_step = True
                    deprotection_depth = max(deprotection_depth, depth)
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Found Boc deprotection step at depth {depth}: {rsmi}")
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if:
    # 1. We have a Boc-protected nitrogen heterocycle intermediate
    # 2. We have both protection and deprotection steps or at least a deprotection step

    # In retrosynthesis, protection should be at higher depth than deprotection
    # (protection happens earlier in forward synthesis)
    correct_sequence = (
        protection_depth > deprotection_depth
        if (has_protection_step and has_deprotection_step)
        else True
    )

    print(f"Has Boc-protected N-heterocycle: {has_boc_protected_n_heterocycle}")
    print(f"Has protection step: {has_protection_step} (depth: {protection_depth})")
    print(f"Has deprotection step: {has_deprotection_step} (depth: {deprotection_depth})")
    print(f"Correct sequence: {correct_sequence}")

    # We need either a protection step or a protected intermediate, plus a deprotection step
    strategy_present = (
        has_boc_protected_n_heterocycle or has_protection_step
    ) and has_deprotection_step

    # If we have both protection and deprotection, ensure correct sequence
    if has_protection_step and has_deprotection_step:
        if correct_sequence:
            if {"type": "sequence", "details": {"before": "Boc protection reaction", "after": "Boc deprotection reaction"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "Boc protection reaction", "after": "Boc deprotection reaction"}})
        strategy_present = strategy_present and correct_sequence

    if has_deprotection_step and has_boc_protected_n_heterocycle:
        if {"type": "co-occurrence", "details": {"targets": ["Boc deprotection reaction", "Boc-protected nitrogen heterocycle"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Boc deprotection reaction", "Boc-protected nitrogen heterocycle"]}})

    return strategy_present, findings_json
