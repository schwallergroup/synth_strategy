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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves a two-step sequence where an α-bromoketone is formed as an intermediate
    and is subsequently consumed in a later step to form a thiazole ring.
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

    has_bromoketone_intermediate = False
    bromoketone_used_for_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal has_bromoketone_intermediate, bromoketone_used_for_heterocycle, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and all(reactants):
                    bromoketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-C-[Br]")
                    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                    # Check if this reaction produces an α-bromoketone
                    if product.HasSubstructMatch(bromoketone_pattern) and not any(
                        r.HasSubstructMatch(bromoketone_pattern) for r in reactants
                    ):
                        print(f"Detected α-bromoketone formation at depth {depth}")
                        has_bromoketone_intermediate = True
                        if "alpha-bromoketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("alpha-bromoketone")
                        if "alpha-bromoketone_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("alpha-bromoketone_formation")

                    # Check if this reaction uses an α-bromoketone to form a thiazole ring
                    if any(
                        r.HasSubstructMatch(bromoketone_pattern) for r in reactants
                    ) and product.HasSubstructMatch(thiazole_pattern) and not any(r.HasSubstructMatch(thiazole_pattern) for r in reactants):
                        print(
                            f"Detected α-bromoketone used for heterocycle formation at depth {depth}"
                        )
                        bromoketone_used_for_heterocycle = True
                        if "alpha-bromoketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("alpha-bromoketone")
                        if "thiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("thiazole")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children with modified depth calculation
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is detected if both conditions are met
    strategy_detected = has_bromoketone_intermediate and bromoketone_used_for_heterocycle
    if strategy_detected:
        print("Complete α-bromoketone intermediate strategy detected")
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "alpha-bromoketone_formation",
                    "thiazole_ring_formation_from_alpha-bromoketone"
                ]
            }
        })

    return strategy_detected, findings_json
