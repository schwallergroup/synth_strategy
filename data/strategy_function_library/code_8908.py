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
    Detects if the synthesis involves a late-stage nucleophilic aromatic substitution
    on a pyrimidine ring (depth 0 or 1).
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

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if one reactant contains a halopyrimidine
            pyrimidine_pattern = Chem.MolFromSmarts("[n]1[c][n][c]([Cl,F,Br,I])[n][c]1")
            aryl_halide_pattern = Chem.MolFromSmarts("c[Cl,F,Br,I]")

            has_pyrimidine = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(pyrimidine_pattern):
                            has_pyrimidine = True
                            if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            if "aryl_halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("aryl_halide")
                except:
                    continue

            # Check if product has a C-N bond where the halogen was
            if has_pyrimidine:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[n]1[c][n][c]([N])[n][c]1")
                    ):
                        result = True
                        if "nucleophilic_aromatic_substitution" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("nucleophilic_aromatic_substitution")
                        if "aryl_amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("aryl_amine")
                        # Add structural constraint if all conditions met
                        if {"type": "positional", "details": {"target": "nucleophilic_aromatic_substitution", "position": "late_stage", "max_depth": 1}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "nucleophilic_aromatic_substitution", "position": "late_stage", "max_depth": 1}})
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result, findings_json