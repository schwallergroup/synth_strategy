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


# Refactoring for Enumeration: Isolate lists into module-level constants
NITROGEN_HETEROCYCLES_OF_INTEREST = [
    "quinoline",
    "carbazole",
    "indole",
    "isoquinoline",
    "pyridine",
    "pyrrole",
]

HETEROCYCLE_PRIORITY_ORDER = [
    "carbazole",
    "isoquinoline",
    "quinoline",
    "indole",
    "pyridine",
    "pyrrole",
]

# Corrected list of relevant named reactions
HETEROCYCLE_FORMATION_REACTIONS = [
    "Fischer indole",
    "Friedlaender chinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation or interconversion of specific nitrogen heterocycles. This is identified in two ways: 1) A reaction step where the set of heterocycles in the product(s) is different from the set in the reactant(s). 2) The reaction is a known named reaction that forms a heterocycle. The specific heterocycles are defined in `NITROGEN_HETEROCYCLES_OF_INTEREST` and the named reactions in `HETEROCYCLE_FORMATION_REACTIONS`.
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

    # Track reaction nodes
    connecting_reactions = []

    def dfs_traverse(node, depth=0):
        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            connecting_reactions.append(node)

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check for direct transformations or known formations in each reaction
    for reaction in connecting_reactions:
        try:
            rsmi = reaction["metadata"]["rsmi"]

            # Check for specific named reactions that form heterocycles
            for reaction_name in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    result = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    # If a named reaction is found, we can return early with the finding
                    return result, findings_json

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check which heterocycles are in reactants (prioritizing larger structures)
            reactant_heterocycles = set()
            for reactant in reactants:
                found_heterocycles = []
                for ring_type in NITROGEN_HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring_type, reactant):
                        found_heterocycles.append(ring_type)
                if found_heterocycles:
                    for ring in HETEROCYCLE_PRIORITY_ORDER:
                        if ring in found_heterocycles:
                            reactant_heterocycles.add(ring)
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                            break

            # Check which heterocycles are in product (prioritizing larger structures)
            product_heterocycles = set()
            found_heterocycles = []
            for ring_type in NITROGEN_HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring_type, product):
                    found_heterocycles.append(ring_type)
            if found_heterocycles:
                for ring in HETEROCYCLE_PRIORITY_ORDER:
                    if ring in found_heterocycles:
                        product_heterocycles.add(ring)
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                        break

            # A transformation occurs if a new heterocycle appears in the product
            if (product_heterocycles and not product_heterocycles.issubset(reactant_heterocycles)):
                result = True
                # No specific structural constraint to add here based on the provided JSON
                # The finding is implicitly captured by the presence of new ring systems in atomic_checks
                return result, findings_json

        except Exception:
            # Silently ignore errors in reaction processing
            continue

    return result, findings_json