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
    "pyrazole", "pyridine", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrimidine", "pyrazine", "triazole", "tetrazole", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole", "piperidine", "pyrrolidine",
    "morpholine", "thiomorpholine", "piperazine", "pyrrole", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole", "quinoline",
    "isoquinoline", "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if any molecule in the synthetic route contains at least two trifluoromethyl groups
    and at least one of the specified heterocyclic scaffolds.
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

    # Track if we found a molecule with multiple fluorinated groups on heterocyclic scaffolds
    strategy_found = False

    def dfs_traverse(node, depth=0):
        nonlocal strategy_found, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Skip if we've already found a match
            if strategy_found:
                return

            # Count trifluoromethyl groups
            trifluoro_indices = (
                checker.get_fg_atom_indices("Trifluoro group", mol_smiles) or []
            )
            fluorinated_groups = len(trifluoro_indices)

            if fluorinated_groups > 0:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # Check for heterocycles
            has_heterocycle = False
            found_heterocycles = []
            if fluorinated_groups >= 2:
                # Record the count constraint if met
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "Trifluoro group",
                        "operator": ">=",
                        "value": 2
                    }
                })

                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, mol_smiles):
                        has_heterocycle = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        found_heterocycles.append(ring)
                        # No break here, as we want to record all found heterocycles

            # Check if this molecule meets our criteria
            if fluorinated_groups >= 2 and has_heterocycle:
                strategy_found = True
                # Record the co-occurrence constraint if met
                co_occurrence_targets = ["Trifluoro group"]
                co_occurrence_targets.extend(found_heterocycles)
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "comment": "The strategy requires the co-occurrence of a Trifluoro group and any one of the ring systems listed in atomic_checks within the same molecule. The list of targets represents this AND/OR logic.",
                        "targets": co_occurrence_targets
                    }
                })

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' (or 'mol'), depth increases
            next_depth = depth + 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return strategy_found, findings_json
