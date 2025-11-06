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


# Refactoring for Enumeration: Isolate the list of heterocycles
HETEROCYCLE_RINGS_OF_INTEREST = [
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
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where nitrile-containing intermediates are used to construct specific heterocyclic scaffolds in early-to-mid stage reactions. The heterocycles checked for are defined in the HETEROCYCLE_RINGS_OF_INTEREST list.
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

    nitrile_intermediates = False
    heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_intermediates, heterocycle_formation, findings_json

        node["depth"] = depth

        if node["type"] == "mol":
            # Check for nitrile groups in intermediates (early stages)
            if depth >= 2:  # Early stages (any step before the final one)
                if checker.check_fg("Nitrile", node["smiles"]):
                    nitrile_intermediates = True
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    # Record the positional constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Nitrile",
                            "condition": "depth >= 2"
                        }
                    })

        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_has_nitrile = False
            for r in reactants:
                if checker.check_fg("Nitrile", r):
                    reactant_has_nitrile = True
                    # No need to add Nitrile to atomic_checks here again, it's covered by mol check
                    break

            product_has_heterocycle = False
            for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                if checker.check_ring(ring, product):
                    product_has_heterocycle = True
                    findings_json["atomic_checks"]["ring_systems"].append(ring)

            reactants_have_heterocycle = False
            for r in reactants:
                for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                    if checker.check_ring(ring, r):
                        reactants_have_heterocycle = True
                        break
                if reactants_have_heterocycle:
                    break

            # Check if this is a heterocycle formation involving a nitrile reactant
            if (
                product_has_heterocycle
                and not reactants_have_heterocycle
                and reactant_has_nitrile
            ):
                heterocycle_formation = True
                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                # Record the co-occurrence constraint if met
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "ring_formation",
                            "Nitrile"
                        ],
                        "scope": "reaction_with_reactant"
                    }
                })

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions were met
    result = nitrile_intermediates and heterocycle_formation
    return result, findings_json
