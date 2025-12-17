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

HETEROCYCLIC_ACID_RINGS = [
    "furan",
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
    "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects an amide coupling reaction where one of the reactants is a carboxylic acid
    bearing a specific heterocyclic ring from a predefined list.
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

    # Track if we've found the required components in a single reaction
    found_valid_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_valid_reaction, findings_json

        if node["type"] == "reaction" and not found_valid_reaction:
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for amide coupling reaction
            is_amide_coupling = False
            amide_coupling_reactions = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Carboxylic acid with primary amine to amide",
                "Schotten-Baumann_amide"
            ]
            for reaction_name in amide_coupling_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    is_amide_coupling = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            if is_amide_coupling:
                # Check for heterocyclic acid
                heterocyclic_acid_found = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", reactant):
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                        # Check if it contains a heterocyclic ring
                        for ring in HETEROCYCLIC_ACID_RINGS:
                            if checker.check_ring(ring, reactant):
                                heterocyclic_acid_found = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break
                    if heterocyclic_acid_found:
                        break

                # If the condition is met in the same reaction
                if heterocyclic_acid_found:
                    found_valid_reaction = True
                    # Add the structural constraint if all conditions are met
                    structural_constraint_obj = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Carboxylic acid with primary amine to amide",
                                "Schotten-Baumann_amide",
                                "Carboxylic acid",
                                "furan",
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
                                "benzimidazole"
                            ]
                        }
                    }
                    if structural_constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraint_obj)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    return found_valid_reaction, findings_json