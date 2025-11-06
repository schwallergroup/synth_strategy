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
    Detects if the synthesis involves a nucleophilic aromatic substitution on a pyrimidine ring.
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

    nas_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nas_found, findings_json

        if node["type"] == "reaction":  # Late in synthesis (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if any reactant has a pyrimidine ring with a halide
                    reactant_has_pyrimidine_halide = False
                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant):
                            if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                            if checker.check_fg("Aromatic halide", reactant):
                                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                                print(f"Found pyrimidine with aromatic halide in reactant: {reactant}")
                                reactant_has_pyrimidine_halide = True
                                break

                    # Check if product has a pyrimidine ring
                    product_has_pyrimidine = checker.check_ring("pyrimidine", product)
                    if product_has_pyrimidine:
                        if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")

                    # Check if this is a nucleophilic aromatic substitution reaction
                    is_nas_reaction = False
                    detected_nas_reaction_name = None

                    nas_reactions = [
                        {"name": "nucl_sub_aromatic_ortho_nitro", "constraint": {"type": "co-occurrence", "details": {"targets": ["pyrimidine", "Aromatic halide", "nucl_sub_aromatic_ortho_nitro"]}}},
                        {"name": "nucl_sub_aromatic_para_nitro", "constraint": {"type": "co-occurrence", "details": {"targets": ["pyrimidine", "Aromatic halide", "nucl_sub_aromatic_para_nitro"]}}},
                        {"name": "heteroaromatic_nuc_sub", "constraint": {"type": "co-occurrence", "details": {"targets": ["pyrimidine", "Aromatic halide", "heteroaromatic_nuc_sub"]}}}
                    ]

                    for reaction_info in nas_reactions:
                        if checker.check_reaction(reaction_info["name"], rsmi):
                            if reaction_info["name"] not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_info["name"])
                            is_nas_reaction = True
                            detected_nas_reaction_name = reaction_info["name"]
                            break

                    if (
                        reactant_has_pyrimidine_halide
                        and product_has_pyrimidine
                        and is_nas_reaction
                    ):
                        print(
                            f"Nucleophilic aromatic substitution on pyrimidine detected at depth {depth}"
                        )
                        nas_found = True
                        # Add the specific structural constraint that was met
                        for reaction_info in nas_reactions:
                            if reaction_info["name"] == detected_nas_reaction_name:
                                if reaction_info["constraint"] not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(reaction_info["constraint"])
                                break
                    else:
                        if not reactant_has_pyrimidine_halide:
                            print(f"No pyrimidine with halide found in reactants at depth {depth}")
                        if not product_has_pyrimidine:
                            print(f"No pyrimidine found in product at depth {depth}")
                        if not is_nas_reaction:
                            print(
                                f"Reaction at depth {depth} is not a nucleophilic aromatic substitution"
                            )

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return nas_found, findings_json
