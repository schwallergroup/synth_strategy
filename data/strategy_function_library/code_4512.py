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


# List of heterocycles that can be formed from 1,3-dicarbonyls
HETEROCYCLES_FROM_DICARBONYLS = [
    "pyrazole",
    "triazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "isoxazole",
    "isothiazole",
    "pyrimidine",
    "pyridazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step strategy where a Claisen-type condensation first forms a 1,3-dicarbonyl compound,
    which is subsequently used to synthesize a specific heterocycle. The target heterocycles are defined
    in the HETEROCYCLES_FROM_DICARBONYLS list.
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

    # Track key transformations
    claisen_condensation = False
    heterocycle_formation = False

    # Track depths
    claisen_depth = -1
    heterocycle_depth = -1

    # Define 1,3-dicarbonyl pattern for specific arrangement check
    dicarbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX3](=[OX1])")

    def dfs_traverse(node, depth=0):
        nonlocal claisen_condensation, heterocycle_formation
        nonlocal claisen_depth, heterocycle_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Check for Claisen condensation (formation of 1,3-dicarbonyl)
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(dicarbonyl_pattern):
                # Check if reactants are esters, ketones, or aldehydes
                ester_count = 0
                ketone_count = 0
                aldehyde_count = 0

                for r in reactants:
                    if checker.check_fg("Ester", r):
                        ester_count += 1
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                    if checker.check_fg("Ketone", r):
                        ketone_count += 1
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if checker.check_fg("Aldehyde", r):
                        aldehyde_count += 1
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

                # Claisen typically needs either two esters or an ester and a carbonyl compound
                if (ester_count >= 2) or (
                    ester_count >= 1 and (ketone_count >= 1 or aldehyde_count >= 1)
                ):
                    claisen_condensation = True
                    claisen_depth = depth
                    if "Claisen condensation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Claisen condensation")

            # Check for heterocycle formation from 1,3-dicarbonyl
            found_heterocycle = False
            for ring in HETEROCYCLES_FROM_DICARBONYLS:
                if checker.check_ring(ring, product):
                    found_heterocycle = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break

            if found_heterocycle:
                # Check if any reactant contains a 1,3-dicarbonyl pattern
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(dicarbonyl_pattern):
                        heterocycle_formation = True
                        heterocycle_depth = depth
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        break

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    result = False
    # Check if we found both transformations in the correct order
    if claisen_condensation and heterocycle_formation:
        # Verify Claisen happens before heterocycle formation (higher depth in retrosynthesis)
        if claisen_depth > heterocycle_depth:
            result = True
            # Add structural constraints if both conditions are met and in order
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Claisen condensation",
                        "ring_formation"
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "Claisen condensation",
                        "ring_formation"
                    ]
                }
            })

    return result, findings_json
