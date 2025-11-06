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


HETEROCYCLE_PATTERNS_OF_INTEREST = {
    "indole": Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]2:[#7]:[#6]:[#6]:[#6]:[#6]2:[#6]1"),
    "pyrazole": Chem.MolFromSmarts("[#7]1:[#7]:[#6]:[#6]:[#6]1"),
    "thiadiazole": Chem.MolFromSmarts("[#7]1:[#7]:[#6]:[#16]:[#6]1"),
    "pyrimidine": Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]:[#6]:[#6]1"),
    "thiazole": Chem.MolFromSmarts("[#6]1:[#6]:[#7]:[#6]:[#16]1"),
    "imidazole": Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]:[#6]1"),
}

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy involving multiple heterocyclic systems
    (at least 3 distinct heterocycles).
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

    # Track heterocycles found in the synthesis
    heterocycles_found = set()

    # Track if the synthesis is convergent
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, heterocycles_found, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                for name, pattern in HETEROCYCLE_PATTERNS_OF_INTEREST.items():
                    if mol.HasSubstructMatch(pattern) and name not in heterocycles_found:
                        heterocycles_found.add(name)
                        findings_json["atomic_checks"]["ring_systems"].append(name)

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if we have multiple reactants (convergent step)
            if len(reactants_smiles) >= 2:
                is_convergent = True
                if "convergent_step" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("convergent_step")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 3 distinct heterocycles and a convergent synthesis
    result = len(heterocycles_found) >= 3 and is_convergent

    if len(heterocycles_found) >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_ring_systems",
                "operator": ">=",
                "value": 3
            }
        })
    
    if is_convergent:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_step",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json