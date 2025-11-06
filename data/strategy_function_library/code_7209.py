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
HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "furan", "thiophene", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole", "triazole", "tetrazole", "isoxazole", "isothiazole",
    "pyrimidine", "pyrazine", "pyridazine", "quinoline", "isoquinoline",
]

N_ARYLATION_REACTIONS = [
    "Buchwald-Hartwig", "Ullmann-Goldberg", "N-arylation",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation_heterocycles",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis step involves the coupling of at least two heterocyclic
    fragments via a specified N-arylation reaction. The function checks if a
    reaction is one of the types listed in `N_ARYLATION_REACTIONS` and if at
    least two of the reactants contain a ring from the
    `HETEROCYCLIC_RINGS_OF_INTEREST` list.
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

    heterocycle_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_coupling_found, findings_json
        if heterocycle_coupling_found: # Optimization: exit early if found
            return

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) >= 2:
                    # Check if the reaction is a known N-arylation type.
                    is_n_arylation = False
                    detected_n_arylations = []
                    for rxn in N_ARYLATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_n_arylation = True
                            detected_n_arylations.append(rxn)
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    if is_n_arylation:
                        # Identify reactants that contain one of the specified heterocycles.
                        heterocycle_reactants = []
                        for r in reactants:
                            detected_rings_for_reactant = []
                            for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                                if checker.check_ring(ring, r):
                                    detected_rings_for_reactant.append(ring)
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                            if detected_rings_for_reactant:
                                heterocycle_reactants.append(r)

                        # The strategy is confirmed if an N-arylation reaction
                        # consumes at least two heterocyclic fragments.
                        if len(heterocycle_reactants) >= 2:
                            heterocycle_coupling_found = True
                            # Add the structural constraint if the condition is met
                            structural_constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "N-arylation reaction",
                                        "coupling of >= 2 heterocyclic reactants"
                                    ],
                                    "scope": "single_reaction_step"
                                }
                            }
                            if structural_constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_obj)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return heterocycle_coupling_found, findings_json
