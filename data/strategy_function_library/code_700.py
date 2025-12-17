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
    This function detects an early-stage aryl-nitrogen coupling strategy,
    where an aryl group is coupled with a nitrogen heterocycle early in the synthesis.
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

    found_early_aryl_n_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_early_aryl_n_coupling, findings_json

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth >= 2)
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl-nitrogen coupling
                aryl_f_pattern = Chem.MolFromSmarts("[c][F]")
                aryl_n_aryl_pattern = Chem.MolFromSmarts("[c][n][c]")

                has_aryl_f = False
                for r in reactants:
                    if r and Chem.MolFromSmiles(r).HasSubstructMatch(aryl_f_pattern):
                        has_aryl_f = True
                        findings_json["atomic_checks"]["functional_groups"].append("aryl fluoride")
                        break

                has_aryl_n_aryl_product = False
                if product:
                    if Chem.MolFromSmiles(product).HasSubstructMatch(aryl_n_aryl_pattern):
                        has_aryl_n_aryl_product = True
                        findings_json["atomic_checks"]["functional_groups"].append("di-aryl amine linkage")

                if has_aryl_f and has_aryl_n_aryl_product:
                    print(f"Found early-stage aryl-nitrogen coupling at depth {depth}")
                    found_early_aryl_n_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append("aryl-nitrogen coupling")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "aryl-nitrogen coupling",
                            "position": "early_stage"
                        }
                    })

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_early_aryl_n_coupling, findings_json