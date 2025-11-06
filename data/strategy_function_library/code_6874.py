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


# Refactoring for Enumeration: Isolate the lists of chemical entities.
NITROGEN_NUCLEOPHILE_RINGS = [
    "piperazine", "pyrrolidine", "morpholine", "piperidine",
    "azetidine", "azepane", "diazepane"
]

AROMATIC_HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "furan", "thiophene",
    "pyrrole", "imidazole", "oxazole", "thiazole", "triazole", "tetrazole",
    "isoxazole", "isothiazole", "oxadiazole", "thiadiazole", "quinoline",
    "isoquinoline", "indole", "benzimidazole", "benzoxazole", "benzothiazole"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage Nucleophilic Aromatic Substitution (SNAr) where a nitrogen nucleophile
    displaces a halogen on an aromatic heterocycle.
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

    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected, findings_json

        if node["type"] == "reaction" and depth <= 2:
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            if len(reactants_smiles) < 2:
                return

            # Check if this is an SNAr reaction using predefined reaction types.
            # Metal-catalyzed couplings (Buchwald, Ullmann) are removed to focus on true SNAr.
            is_snar = False
            snar_reaction_types = [
                "heteroaromatic_nuc_sub",
                "nucl_sub_aromatic_ortho_nitro",
                "nucl_sub_aromatic_para_nitro",
                "N-arylation_heterocycles"
            ]
            for r_name in snar_reaction_types:
                if checker.check_reaction(r_name, rsmi):
                    is_snar = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)

            # If a general SNAr-type reaction is detected, verify it involves the specific
            # components of interest: a nitrogen nucleophile and a halogenated heterocycle.
            # The original manual check fallback was removed due to a critical flaw in its
            # product verification logic, which would cause false positives/negatives.
            if is_snar:
                nitrogen_nucleophile_found = False
                heterocycle_with_halogen_found = False

                for reactant in reactants_smiles:
                    # Check for primary/secondary amine nucleophiles, including specified cyclic amines.
                    # Tertiary amines are correctly excluded as they are not nucleophilic.
                    if checker.check_fg("Primary amine", reactant):
                        nitrogen_nucleophile_found = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", reactant):
                        nitrogen_nucleophile_found = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                    for ring in NITROGEN_NUCLEOPHILE_RINGS:
                        if checker.check_ring(ring, reactant):
                            nitrogen_nucleophile_found = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    # Check for a halogen on one of the specified aromatic heterocycles.
                    if checker.check_fg("Aromatic halide", reactant):
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        for ring_type in AROMATIC_HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring_type, reactant):
                                heterocycle_with_halogen_found = True
                                if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                
                # If both specific components are present in a confirmed SNAr-type reaction, flag it.
                if nitrogen_nucleophile_found and heterocycle_with_halogen_found:
                    snar_detected = True
                    # Add structural constraint for co-occurrence
                    if {"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["SNAr_reaction_type", "nitrogen_nucleophile_reactant", "halogenated_heterocycle_reactant"]}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["SNAr_reaction_type", "nitrogen_nucleophile_reactant", "halogenated_heterocycle_reactant"]}})
                    # Add structural constraint for positional (late-stage)
                    if depth <= 2:
                        if {"type": "positional", "details": {"target": "SNAr reaction", "position": "within last 3 stages"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "SNAr reaction", "position": "within last 3 stages"}})
                    return

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return snar_detected, findings_json
