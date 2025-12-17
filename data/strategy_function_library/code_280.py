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
    "pyridine", "pyrrole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrazole", "pyrimidine", "piperidine", "morpholine",
    "piperazine", "indole", "benzimidazole", "quinoline", "isoquinoline",
    "pyrazine", "pyridazine"
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving late-stage amide coupling
    between two fragments, each containing a specific heterocycle.
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

    found_amide_coupling = False
    is_final_step = False

    first_reaction = None
    if route["type"] == "mol" and route.get("children"):
        for child in route.get("children", []):
            if child["type"] == "reaction":
                first_reaction = child
                break

    def dfs_traverse(node, current_depth=0):
        nonlocal found_amide_coupling, is_final_step, findings_json

        if node["type"] == "reaction":
            depth = node.get("metadata", {}).get("depth", current_depth)

            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi and ">" in rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_amide_coupling = False

                if product_smiles and len(reactants_smiles) >= 2:
                    acid_reactant = None
                    for r in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r):
                            acid_reactant = r
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            break

                    amine_reactant = None
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amine", r):
                            amine_reactant = r
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            break
                        if checker.check_fg("Secondary amine", r):
                            amine_reactant = r
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            break

                    has_amide = False
                    if checker.check_fg("Primary amide", product_smiles):
                        has_amide = True
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if checker.check_fg("Secondary amide", product_smiles):
                        has_amide = True
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product_smiles):
                        has_amide = True
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    heterocyclic_reactants = 0
                    for r in reactants_smiles:
                        for ring in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring, r):
                                heterocyclic_reactants += 1
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break

                    if (
                        acid_reactant
                        and amine_reactant
                        and has_amide
                        and heterocyclic_reactants >= 2
                    ):
                        is_amide_coupling = True
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                        if {"type": "count", "details": {"target": "reactant_with_specified_heterocycle", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactant_with_specified_heterocycle", "operator": ">=", "value": 2}})

                if is_amide_coupling:
                    found_amide_coupling = True
                    if depth == 1:
                        is_final_step = True
                        if {"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}})

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, current_depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    result = found_amide_coupling and is_final_step
    return result, findings_json
