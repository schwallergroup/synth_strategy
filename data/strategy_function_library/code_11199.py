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


ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid to amide conversion",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "{Schotten-Baumann_amide}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step synthetic sequence: (1) ester hydrolysis to a carboxylic acid, followed by (2) amide formation from the acid. The specific reaction types checked for each step are defined in the module-level constants ESTER_HYDROLYSIS_REACTIONS and AMIDE_FORMATION_REACTIONS.
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

    transformations_by_depth = {}
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester to acid transformation
                ester_found_in_reactant = False
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant):
                        ester_found_in_reactant = True
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")

                carboxylic_acid_found_in_product = checker.check_fg("Carboxylic acid", product)
                if carboxylic_acid_found_in_product:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                if ester_found_in_reactant and carboxylic_acid_found_in_product:
                    ester_hydrolysis_reaction_found = False
                    for reaction_type in ESTER_HYDROLYSIS_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            ester_hydrolysis_reaction_found = True
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    if ester_hydrolysis_reaction_found:
                        transformations_by_depth[depth] = {
                            "type": "ester_to_acid",
                        }

                # Check for acid to amide transformation
                carboxylic_acid_found_in_reactant = False
                for reactant in reactants:
                    if checker.check_fg("Carboxylic acid", reactant):
                        carboxylic_acid_found_in_reactant = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                amide_found_in_product = False
                if checker.check_fg("Primary amide", product):
                    amide_found_in_product = True
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product):
                    amide_found_in_product = True
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product):
                    amide_found_in_product = True
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                if carboxylic_acid_found_in_reactant and amide_found_in_product:
                    amide_formation_reaction_found = False
                    for reaction_type in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            amide_formation_reaction_found = True
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    if amide_formation_reaction_found:
                        transformations_by_depth[depth] = {
                            "type": "acid_to_amide",
                        }

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check for direct transformations in the correct sequence.
    # An acid_to_amide step must occur at a lower depth (later in synthesis)
    # than an ester_to_acid step (d1 < d2).
    for d1 in transformations_by_depth:
        if transformations_by_depth[d1].get("type") == "acid_to_amide":
            for d2 in transformations_by_depth:
                if d2 > d1 and transformations_by_depth[d2].get("type") == "ester_to_acid":
                    result = True
                    # Add the structural constraint if found
                    structural_constraint_obj = {
                        "type": "sequence",
                        "details": {
                            "ordered_events": [
                                "ester_hydrolysis",
                                "amide_formation"
                            ],
                            "description": "An ester hydrolysis reaction must occur at an earlier synthetic step (greater depth) than an amide formation reaction."
                        }
                    }
                    if structural_constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraint_obj)

    return result, findings_json
