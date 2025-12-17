from typing import Tuple, Dict, List
import copy
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
    """Detects the formation of a pyrazole ring in the final synthetic step (depth=1) from a hydrazine or hydrazone precursor. The strategy is confirmed if a pyrazole ring is present in the product but not in any reactant, and at least one reactant contains a hydrazine or hydrazone functional group."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    pyrazole_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_formation_found, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if depth == 1:
                    product_has_pyrazole = checker.check_ring("pyrazole", product_smiles)
                    reactant_has_pyrazole = any(
                        checker.check_ring("pyrazole", r) for r in reactants_smiles
                    )

                    if product_has_pyrazole:
                        if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

                    if product_has_pyrazole and not reactant_has_pyrazole:
                        hydrazine_found = False
                        for r in reactants_smiles:
                            if checker.check_fg("Hydrazine", r):
                                if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                                hydrazine_found = True
                            if checker.check_fg("Hydrazone", r):
                                if "Hydrazone" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Hydrazone")
                                hydrazine_found = True

                        if hydrazine_found:
                            pyrazole_formation_found = True
                            # Add structural constraints when the main condition is met
                            if {"type": "positional", "details": {"target": "ring_formation", "event_type": "named_reactions", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "event_type": "named_reactions", "position": "last_stage"}})
                            if {"type": "co-occurrence", "details": {"scope": "single_reaction", "targets": ["ring_formation", ["Hydrazine", "Hydrazone"]]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "single_reaction", "targets": ["ring_formation", ["Hydrazine", "Hydrazone"]]}})

            except Exception:
                # Ignore errors in reaction analysis to ensure robustness
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return pyrazole_formation_found, findings_json