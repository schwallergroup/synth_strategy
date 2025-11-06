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
HETEROCYCLES_OF_INTEREST = [
    "furan", "pyrrole", "thiophene", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine",
    "triazole", "tetrazole", "indole", "benzofuran", "benzothiophene",
    "quinoline", "isoquinoline", "benzimidazole", "benzoxazole", "benzothiazole",
]

ALKYNE_COUPLING_REACTIONS = [
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf", "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_alkenyl halide", "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl OTf", "Sonogashira alkyne_alkenyl OTf",
    "Sonogashira acetylene_acyl halide", "Sonogashira alkyne_acyl halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage Sonogashira coupling that forms an alkyne linker
    between two different heterocyclic systems. The specific heterocycles and
    Sonogashira variants are enumerated in the module-level lists
    HETEROCYCLES_OF_INTEREST and ALKYNE_COUPLING_REACTIONS.
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

    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, findings_json

        # Check for late-stage coupling (depth 0-2)
        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactant_smiles = reactants_part.split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                is_coupling = False
                for reaction in ALKYNE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction, rsmi):
                        is_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction)
                        break

                # Check for alkyne in product
                has_alkyne_product = checker.check_fg("Alkyne", product_smiles)
                if has_alkyne_product:
                    findings_json["atomic_checks"]["functional_groups"].append("Alkyne")

                # Check for heterocycles in reactants
                reactant_heterocycles = []
                for i, reactant in enumerate(reactant_smiles):
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.append((i, heterocycle))
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break # Found a heterocycle, move to next reactant

                # Check for heterocycles in product
                product_heterocycles = []
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.append(heterocycle)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # Verify we have the pattern:
                # 1. It's a specified coupling reaction
                # 2. Product has an alkyne
                # 3. Product has at least 2 heterocycles
                # 4. Reactants had heterocycles in at least 2 different molecules
                if (
                    is_coupling
                    and has_alkyne_product
                    and len(product_heterocycles) >= 2
                    and len(set([i for i, _ in reactant_heterocycles])) >= 2
                ):
                    found_pattern = True
                    # Record structural constraints
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Sonogashira coupling",
                            "constraint": {
                                "on": "depth",
                                "operator": "<=",
                                "value": 2
                            }
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Sonogashira coupling",
                                "product_with_alkyne"
                            ],
                            "scope": "single_reaction_step"
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "specified_heterocycles_in_product",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactant_molecules_with_specified_heterocycle",
                            "operator": ">=",
                            "value": 2
                        }
                    })
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            if not found_pattern:  # Stop traversal if pattern already found
                # New logic for depth calculation
                if node["type"] == "reaction":
                    # If current node is a reaction, depth remains the same for children
                    dfs_traverse(child, depth)
                else:
                    # If current node is chemical, depth increases for children
                    dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_pattern, findings_json
