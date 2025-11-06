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


HETEROCYCLES_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran", "oxirane",
    "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene", "pyrrole", "pyridine",
    "pyrazole", "imidazole", "oxazole", "thiazole", "pyrimidine", "pyrazine",
    "pyridazine", "triazole", "tetrazole", "indole", "quinoline", "isoquinoline",
    "benzoxazole", "benzothiazole", "benzimidazole", "thiophene", "thiopyran",
    "thiirane", "thietane",
]

def main(route) -> Tuple[bool, Dict]:
    """Checks for the de novo formation of specific heterocyclic rings, as defined in the HETEROCYCLES_OF_INTEREST list (e.g., furan, pyridine, indole). A ring is considered formed if it is present in the product but absent from all reactants."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    heterocyclic_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocyclic_formation_found, findings_json
        if heterocyclic_formation_found:  # Early exit to avoid unnecessary computation
            return

        if node.get("type") == "reaction" and "rsmi" in node.get("metadata", {}):
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part, product = rsmi.split(">>", 1)
                reactants = reactants_part.split(".")

                for ring in HETEROCYCLES_OF_INTEREST:
                    # This check is robust: it identifies if a specific heterocyclic ring system
                    # is created in this step (i.e., present in product, absent in all reactants).
                    # This avoids false positives from simple group transfers and is less ambiguous
                    # than attempting to count rings with a boolean checker.
                    if checker.check_ring(ring, product) and not any(checker.check_ring(ring, r) for r in reactants):
                        heterocyclic_formation_found = True
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                        # Add the structural constraint if a ring formation is found
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "ring_formation",
                                "operator": ">=",
                                "value": 1
                            }
                        })
                        return  # Exit after finding one, no need to check other rings

            except (ValueError, KeyError):
                # Fail silently on malformed RSMI or missing metadata
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocyclic_formation_found, findings_json
