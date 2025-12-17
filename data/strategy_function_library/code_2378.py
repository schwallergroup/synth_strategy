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


SNAR_REACTION_TYPES = [
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses where the final product contains both a thiophene and a pyrimidine ring. The strategy is confirmed if at least one reaction step in the synthesis is a nucleophilic aromatic substitution (SNAr) that involves a thiophene-containing reactant and a halopyrimidine reactant. The check for SNAr is performed using a specific list of reaction types: 'heteroaromatic_nuc_sub', 'nucl_sub_aromatic_ortho_nitro', and 'nucl_sub_aromatic_para_nitro'.
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

    # Track if we found SNAr on pyrimidine
    found_snar = False
    # Track if thiophene-pyrimidine system is present in final product
    has_thiophene_pyrimidine = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar, has_thiophene_pyrimidine, findings_json

        # Check if final product has thiophene-pyrimidine system
        if node["type"] == "mol" and depth == 0:
            mol_smiles = node["smiles"]
            thiophene_found = checker.check_ring("thiophene", mol_smiles)
            pyrimidine_found = checker.check_ring("pyrimidine", mol_smiles)

            if thiophene_found:
                if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("thiophene")
            if pyrimidine_found:
                if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")

            if thiophene_found and pyrimidine_found:
                has_thiophene_pyrimidine = True
                # Add structural constraint for final product co-occurrence
                constraint_obj = {
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "thiophene",
                            "pyrimidine"
                        ],
                        "scope": "molecule",
                        "position": "last_stage"
                    }
                }
                if constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint_obj)

        # Check reactions for SNAr
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a nucleophilic aromatic substitution reaction
                is_snar = False
                for rxn in SNAR_REACTION_TYPES:
                    if checker.check_reaction(rxn, rsmi):
                        is_snar = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                if is_snar:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product has pyrimidine and thiophene
                    product_has_pyrimidine = checker.check_ring("pyrimidine", product)
                    product_has_thiophene = checker.check_ring("thiophene", product)

                    if product_has_pyrimidine:
                        if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                    if product_has_thiophene:
                        if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("thiophene")

                    if product_has_pyrimidine and product_has_thiophene:
                        # Check if any reactant has pyrimidine with halide
                        has_halopyrimidine = False
                        for r in reactants:
                            r_has_pyrimidine = checker.check_ring("pyrimidine", r)
                            r_has_aromatic_halide = checker.check_fg("Aromatic halide", r)

                            if r_has_pyrimidine:
                                if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                            if r_has_aromatic_halide:
                                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                            if r_has_pyrimidine and r_has_aromatic_halide:
                                has_halopyrimidine = True
                                break

                        # Check if any reactant has thiophene
                        has_thiophene = False
                        for r in reactants:
                            r_has_thiophene = checker.check_ring("thiophene", r)
                            if r_has_thiophene:
                                has_thiophene = True
                                if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append("thiophene")
                                break

                        if has_halopyrimidine and has_thiophene:
                            found_snar = True
                            # Add structural constraint for reaction co-occurrence
                            constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "thiophene",
                                        "pyrimidine",
                                        "Aromatic halide"
                                    ],
                                    "scope": "reactants_of_specific_reaction",
                                    "reaction_names": SNAR_REACTION_TYPES,
                                    "condition": "The 'pyrimidine' and 'Aromatic halide' must be on the same reactant molecule."
                                }
                            }
                            if constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = found_snar and has_thiophene_pyrimidine
    return result, findings_json
