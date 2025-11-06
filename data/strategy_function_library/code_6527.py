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


AROMATIC_RINGS_OF_INTEREST = [
    "thiophene",
    "benzene",
    "pyridine",
    "furan",
    "pyrrole",
    "indole",
    "naphthalene",
    "quinoline",
    "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage biaryl coupling within the final step (`depth <= 1`). The strategy is identified by finding reactions where: 1) At least two reactants contain an aromatic ring from the `AROMATIC_RINGS_OF_INTEREST` list. 2) At least one reactant contains a Boronic acid, Boronic ester, or Aromatic halide. 3) The product contains an aromatic ring from the `AROMATIC_RINGS_OF_INTEREST` list.
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

    biaryl_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_coupling_detected, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late-stage reaction (depth 0, 1, or 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # The inefficient and redundant named reaction check block has been removed.

                # General approach using a robust chemical heuristic.
                if not biaryl_coupling_detected and depth <= 1:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]
                    product_mol = Chem.MolFromSmiles(product)

                    # Record positional constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "biaryl_coupling_event",
                            "position": "depth <= 1"
                        }
                    })

                    if product_mol:
                        # Check if product contains aromatic rings
                        product_has_aromatic = False
                        for ring in AROMATIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, product):
                                product_has_aromatic = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break

                        if product_has_aromatic:
                            # Count aromatic rings in reactants
                            aromatic_reactants_count = 0
                            for reactant in reactants:
                                for ring in AROMATIC_RINGS_OF_INTEREST:
                                    if checker.check_ring(ring, reactant):
                                        aromatic_reactants_count += 1
                                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                                        break

                            if aromatic_reactants_count >= 2:
                                # Record count constraint if met
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "reactants_with_aromatic_ring_of_interest",
                                        "operator": ">=",
                                        "value": 2,
                                        "scope": "per_reaction"
                                    }
                                })

                                # Additional check: look for boronic acids/esters or halides (common in couplings)
                                boronic_or_halide = False
                                for reactant in reactants:
                                    if checker.check_fg("Boronic acid", reactant):
                                        boronic_or_halide = True
                                        if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                                    if checker.check_fg("Boronic ester", reactant):
                                        boronic_or_halide = True
                                        if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                                    if checker.check_fg("Aromatic halide", reactant):
                                        boronic_or_halide = True
                                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                                    if boronic_or_halide:
                                        break

                                if boronic_or_halide:
                                    print(
                                        f"Late-stage biaryl coupling detected at depth {depth} (general approach)"
                                    )
                                    biaryl_coupling_detected = True
                                    # Record co-occurrence constraint if met
                                    findings_json["structural_constraints"].append({
                                        "type": "co-occurrence",
                                        "details": {
                                            "targets": [
                                                "product_has_aromatic_ring_of_interest",
                                                "reactant_has_coupling_functional_group"
                                            ],
                                            "scope": "per_reaction"
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

    return biaryl_coupling_detected, findings_json
