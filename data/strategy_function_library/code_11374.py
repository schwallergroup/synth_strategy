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


# Refactoring for Enumeration: Isolate the lists of heterocycles
HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazine",
    "furan", "thiophene", "pyrrole", "imidazole", "oxazole", "thiazole",
    "triazole", "tetrazole", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole", "quinoline", "isoquinoline",
]
N_NUCLEOPHILIC_HETEROCYCLES = [
    "pyrrole", "imidazole", "pyrazole", "triazole", "tetrazole", "indole", "benzimidazole"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for a late-stage coupling of two heterocycles via a specific SNAr reaction. The strategy is identified in the final synthetic step if a recognized SNAr reaction occurs between a reactant that is an aromatic halide and contains a ring from the `HETEROCYCLES_OF_INTEREST` list, and a second reactant that is a nucleophilic heterocycle from the `N_NUCLEOPHILIC_HETEROCYCLES` list.
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

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json

        if node["type"] == "reaction" and depth <= 1:
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Check for specific SNAr reaction types
            is_snar = False
            if checker.check_reaction("heteroaromatic_nuc_sub", rsmi):
                is_snar = True
                findings_json["atomic_checks"]["named_reactions"].append("heteroaromatic_nuc_sub")
            elif checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi):
                is_snar = True
                findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_ortho_nitro")
            elif checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi):
                is_snar = True
                findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_para_nitro")

            if not is_snar:
                return

            # Identify the heterocyclic halide and heterocyclic N-nucleophile reactants
            halide_reactant = None
            n_nucleophile_reactant = None

            for reactant in reactants:
                if checker.check_fg("Aromatic halide", reactant):
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                    for ring_name in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring_name, reactant):
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            halide_reactant = reactant
                            break
                if halide_reactant:
                    break

            for reactant in reactants:
                if reactant != halide_reactant:
                    for ring_name in N_NUCLEOPHILIC_HETEROCYCLES:
                        if checker.check_ring(ring_name, reactant):
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            n_nucleophile_reactant = reactant
                            break
                if n_nucleophile_reactant:
                    break

            # Final check: SNAr reaction + both required heterocyclic reactants found
            if halide_reactant and n_nucleophile_reactant:
                result = True
                # Add structural constraints if the main condition is met
                # Assuming 'last_stage' means depth <= 1 for this specific function
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "SNAr_heterocycle_coupling",
                        "position": "last_stage"
                    }
                })
                # This co-occurrence is implicitly checked by the logic that finds both reactants
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "heteroaromatic_nuc_sub",
                            "Aromatic halide"
                        ]
                    }
                })

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return result, findings_json