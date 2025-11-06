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


HETEROCYCLIC_RINGS = [
    "imidazole", "pyrazole", "triazole", "tetrazole", "oxazole",
    "thiazole", "isoxazole", "isothiazole", "pyrrole", "furan",
    "thiophene", "pyridine", "pyrimidine", "pyrazine", "pyridazine",
    "indole", "benzimidazole", "benzoxazole", "benzothiazole",
    "quinoline", "isoquinoline",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Paal-Knorr pyrrole synthesis", "Fischer indole",
    "benzimidazole_derivatives_aldehyde",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzothiazole",
    "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "tetrazole_terminal", "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst", "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester", "pyrazole", "oxadiazole",
    "imidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the *de novo* construction of a heterocyclic scaffold. It first checks if the reaction is a known heterocycle-forming reaction (e.g., Paal-Knorr, Fischer indole) from the HETEROCYCLE_FORMING_REACTIONS list. As a fallback, it verifies the formation of a specific heterocycle from the HETEROCYCLIC_RINGS list in the product, ensuring it was not present in any of the reactants.
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

    has_scaffold_construction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_scaffold_construction, findings_json
        if has_scaffold_construction: # Early exit if already found
            return

        if node["type"] == "reaction" and node.get("reaction_obj"):
            reaction = node["reaction_obj"]

            # Check 1: Use specific named reaction checkers
            for reaction_type in HETEROCYCLE_FORMING_REACTIONS:
                if checker.check_reaction(reaction_type, reaction):
                    has_scaffold_construction = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    return

            # Check 2: Fallback to structural check for ring formation
            product = reaction.product
            reactants = reaction.reactants

            for ring_name in HETEROCYCLIC_RINGS:
                # Check if the ring is formed in the product
                if checker.check_ring(ring_name, product):
                    # And ensure it was not present in any reactant
                    if not any(checker.check_ring(ring_name, r) for r in reactants):
                        has_scaffold_construction = True
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        return

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_scaffold_construction, findings_json
