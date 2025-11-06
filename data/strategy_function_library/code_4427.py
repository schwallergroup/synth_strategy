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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a 2-fluoropyridine moiety is introduced in the late stages of the synthesis.
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

    fluoropyridine_introduction_depth = -1
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal fluoropyridine_introduction_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for fluoropyridine pattern
            fluoropyridine_pattern = Chem.MolFromSmarts("c1ncccc1F")
            pyridine_pattern = Chem.MolFromSmarts("n1ccccc1")

            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return

            # Check if product has fluoropyridine
            if product_mol.HasSubstructMatch(fluoropyridine_pattern):
                if "2-fluoropyridine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("2-fluoropyridine")
                if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                # Check if any reactant has fluoropyridine
                reactant_has_fluoropyridine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(fluoropyridine_pattern):
                        reactant_has_fluoropyridine = True
                        break

                if not reactant_has_fluoropyridine:
                    if (fluoropyridine_introduction_depth == -1 or depth < fluoropyridine_introduction_depth):
                        fluoropyridine_introduction_depth = depth
                        print(f"Fluoropyridine introduction detected at depth {depth}")
                        if "2-fluoropyridine_introduction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("2-fluoropyridine_introduction")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late introduction if depth is less than 3
    is_late_introduction = (
        fluoropyridine_introduction_depth >= 0 and fluoropyridine_introduction_depth < 3
    )
    print(f"Fluoropyridine introduction depth: {fluoropyridine_introduction_depth}")

    result = is_late_introduction

    if result:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "2-fluoropyridine_introduction",
                "position": "late_stage"
            }
        })

    return result, findings_json