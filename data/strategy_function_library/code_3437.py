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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a pyrazole-aldehyde motif is preserved throughout
    the synthesis while other transformations occur around it.
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

    # Track if the final product has pyrazole-aldehyde and if it's preserved
    final_product_has_motif = False
    motif_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_motif, motif_preserved, findings_json

        if node["type"] == "mol" and "smiles" in node and node["smiles"]:
            mol_smiles = node["smiles"]
            has_pyrazole = checker.check_ring("pyrazole", mol_smiles)
            if has_pyrazole:
                if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
            has_aldehyde = checker.check_fg("Aldehyde", mol_smiles)
            if has_aldehyde:
                if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
            has_motif = has_pyrazole and has_aldehyde

            # Check if this is the final product (depth 0)
            if depth == 0:
                final_product_has_motif = has_motif
                print(f"Final product has pyrazole-aldehyde motif: {final_product_has_motif}")

            # Process children (reactions)
            for child in node.get("children", []):
                if (
                    child["type"] == "reaction"
                    and "metadata" in child
                    and "rsmi" in child["metadata"]
                ):
                    rsmi = child["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has both pyrazole and aldehyde
                    reactant_has_motif = False
                    for reactant in reactants:
                        if checker.check_ring("pyrazole", reactant):
                            if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                        if checker.check_fg("Aldehyde", reactant):
                            if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

                        if checker.check_ring("pyrazole", reactant) and checker.check_fg(
                            "Aldehyde", reactant
                        ):
                            reactant_has_motif = True
                            break

                    # If a reactant has the motif but it's modified in the reaction
                    if reactant_has_motif and not (
                        checker.check_ring("pyrazole", product)
                        and checker.check_fg("Aldehyde", product)
                    ):
                        motif_preserved = False
                        # Record the negation constraint if motif is destroyed
                        if {"type": "negation", "details": {"target": "destruction of pyrazole-aldehyde motif"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "destruction of pyrazole-aldehyde motif"}})
                        print(f"Pyrazole-aldehyde motif modified in reaction: {rsmi}")

                # Continue traversal
                # New logic: depth increases only from chemical to reaction
                new_depth = depth + 1 if node["type"] != "reaction" else depth
                dfs_traverse(child, new_depth)
        elif node["type"] == "reaction":
            # Process children of reaction nodes
            for child in node.get("children", []):
                # New logic: depth increases only from chemical to reaction
                new_depth = depth + 1 if node["type"] != "reaction" else depth
                dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # The strategy is successful if the final product has the motif and it's preserved
    result = final_product_has_motif and motif_preserved

    # Record the co-occurrence constraint if the final product has the motif
    if final_product_has_motif:
        co_occurrence_constraint = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "pyrazole",
                    "Aldehyde"
                ],
                "position": "last_stage"
            }
        }
        if co_occurrence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(co_occurrence_constraint)

    return result, findings_json
