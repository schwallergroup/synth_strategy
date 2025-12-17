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


HETEROCYCLES_OF_INTEREST = [
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "furan",
    "thiophene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving the formation of a bromine-containing heterocycle
    that is maintained throughout the synthesis.
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

    found_heterocycle_formation = False
    maintains_bromine = True # This flag is always true in the original logic, so it will remain true.
    heterocycle_with_bromine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation, maintains_bromine, heterocycle_with_bromine_found, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle formation (typically early in synthesis)
                if depth >= 3:  # Early in synthesis (higher depth)
                    reactant_has_heterocycle_with_br = False
                    for reactant in reactants_smiles:
                        if "Br" in reactant:
                            if "bromine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("bromine")
                        for ring_type in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring_type, reactant):
                                if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                                if "Br" in reactant:
                                    reactant_has_heterocycle_with_br = True
                                    break
                        if reactant_has_heterocycle_with_br:
                            break

                    product_has_heterocycle_with_br = False
                    if "Br" in product_smiles:
                        if "bromine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("bromine")
                    for ring_type in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring_type, product_smiles):
                            if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                            if "Br" in product_smiles:
                                product_has_heterocycle_with_br = True
                                break

                    if product_has_heterocycle_with_br and not reactant_has_heterocycle_with_br:
                        found_heterocycle_formation = True
                        # Record the structural constraint
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "bromo_heterocycle_formation",
                                "position": "early_stage",
                                "description": "A bromine-containing heterocycle must be formed (i.e., present in products but not reactants) in an early stage of the synthesis (depth >= 3)."
                            }
                        })

                # Check if any stage has a bromine-containing heterocycle
                if "Br" in product_smiles:
                    if "bromine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("bromine")
                    for ring_type in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring_type, product_smiles):
                            if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                            heterocycle_with_bromine_found = True
                            break
            except Exception as e:
                pass

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = found_heterocycle_formation and maintains_bromine and heterocycle_with_bromine_found
    return result, findings_json
