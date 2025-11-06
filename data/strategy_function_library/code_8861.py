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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route preserves heterocyclic cores
    (benzofuran and thiophene) throughout the synthesis.
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

    # Track heterocycles through the synthesis
    preserved_benzofuran = True
    preserved_thiophene = True
    result = False

    # Check if final product has either heterocycle
    final_product_smiles = route["smiles"]

    # Check for benzofuran
    has_benzofuran_in_product = checker.check_ring("benzofuran", final_product_smiles)
    if has_benzofuran_in_product:
        findings_json["atomic_checks"]["ring_systems"].append("benzofuran")

    # Check for thiophene
    has_thiophene_in_product = checker.check_ring("thiophene", final_product_smiles)
    if has_thiophene_in_product:
        findings_json["atomic_checks"]["ring_systems"].append("thiophene")

    # If neither heterocycle is in the final product, no need to check preservation
    if not has_benzofuran_in_product and not has_thiophene_in_product:
        return False, findings_json

    def dfs_traverse(node, depth=0):
        nonlocal preserved_benzofuran, preserved_thiophene, findings_json

        if node["type"] == "reaction":
            # Get reactants and product from reaction
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if benzofuran is preserved in this reaction
            if has_benzofuran_in_product:
                product_has_benzofuran = checker.check_ring("benzofuran", product_smiles)

                # Check reactants
                reactants_have_benzofuran = False
                for r in reactants_smiles:
                    if checker.check_ring("benzofuran", r):
                        reactants_have_benzofuran = True
                        break

                # If product has benzofuran but reactants don't, it's not preserved
                if product_has_benzofuran and not reactants_have_benzofuran:
                    preserved_benzofuran = False
                    # This implies a ring formation, which is a negation constraint
                    if {"type": "negation", "details": {"target": "ring_formation", "context": "benzofuran"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "ring_formation", "context": "benzofuran"}})

            # Check if thiophene is preserved in this reaction
            if has_thiophene_in_product:
                product_has_thiophene = checker.check_ring("thiophene", product_smiles)

                # Check reactants
                reactants_have_thiophene = False
                for r in reactants_smiles:
                    if checker.check_ring("thiophene", r):
                        reactants_have_thiophene = True
                        break

                # If product has thiophene but reactants don't, it's not preserved
                if product_has_thiophene and not reactants_have_thiophene:
                    preserved_thiophene = False
                    # This implies a ring formation, which is a negation constraint
                    if {"type": "negation", "details": {"target": "ring_formation", "context": "thiophene"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "ring_formation", "context": "thiophene"}})

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return true if at least one heterocycle is in the product and was preserved
    if has_benzofuran_in_product and preserved_benzofuran:
        result = True
        # Add the co-occurrence constraint if benzofuran was found and preserved
        if {"type": "co-occurrence", "details": {"targets": ["benzofuran", "thiophene"], "target_type": "ring_system", "min_occurrences": 1, "scope": "final_product"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["benzofuran", "thiophene"], "target_type": "ring_system", "min_occurrences": 1, "scope": "final_product"}})

    if has_thiophene_in_product and preserved_thiophene:
        result = True
        # Add the co-occurrence constraint if thiophene was found and preserved
        if {"type": "co-occurrence", "details": {"targets": ["benzofuran", "thiophene"], "target_type": "ring_system", "min_occurrences": 1, "scope": "final_product"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["benzofuran", "thiophene"], "target_type": "ring_system", "min_occurrences": 1, "scope": "final_product"}})

    return result, findings_json