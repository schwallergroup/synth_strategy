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


ELECTRON_DEFICIENT_HETEROCYCLES = [
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazine",
    "triazole",
    "tetrazole",
    "oxazole",
    "thiazole",
    "isoxazole",
    "isothiazole",
    "quinoline",
    "isoquinoline",
    "imidazole",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects nucleophilic aromatic substitution (SNAr) on electron-deficient heterocycles,
    specifically halogen displacement by an amine.
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

    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for reactions where an amine and halogenated heterocycle react
            has_heterocycle = False
            has_amine = False
            detected_heterocycle_name = None

            for reactant in reactants:
                # Check for aromatic halide on a specified heterocycle
                if checker.check_fg("Aromatic halide", reactant):
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                    for heterocycle in ELECTRON_DEFICIENT_HETEROCYCLES:
                        if checker.check_ring(heterocycle, reactant):
                            has_heterocycle = True
                            detected_heterocycle_name = heterocycle
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break

                # Check for amine nucleophile
                if checker.check_fg("Primary amine", reactant):
                    has_amine = True
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                elif checker.check_fg("Secondary amine", reactant):
                    has_amine = True
                    if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                elif checker.check_fg("Aniline", reactant):
                    has_amine = True
                    if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aniline")

            # If we have the required components, check the product
            if has_heterocycle and has_amine:
                # The check for a tertiary amine even if a halide is present handles polyhalogenated cases
                # where only one halide is substituted.
                if checker.check_fg("Tertiary amine", product):
                    if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                    found_snar = True
                elif (not checker.check_fg("Aromatic halide", product) and (
                        checker.check_fg("Secondary amine", product) or checker.check_fg("Tertiary amine", product)
                    )):
                    if checker.check_fg("Secondary amine", product):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    if checker.check_fg("Tertiary amine", product):
                        if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                    found_snar = True

                if found_snar:
                    # Add SNAr to named_reactions if not already present
                    if "SNAr" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("SNAr")
                    
                    # Add structural constraint if not already present
                    snar_constraint = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Aromatic halide",
                                "electron_deficient_heterocycle",
                                "amine_nucleophile"
                            ],
                            "scope": "reaction_step",
                            "notes": "The strategy identifies a single reaction step where one reactant has both an 'Aromatic halide' and an 'electron_deficient_heterocycle' (any ring from the ring_systems list), and another reactant is an 'amine_nucleophile' ('Primary amine', 'Secondary amine', or 'Aniline')."
                        }
                    }
                    if snar_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(snar_constraint)

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_snar, findings_json
