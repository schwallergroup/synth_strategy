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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route contains a dimethoxybenzene motif throughout the synthesis.
    This indicates a strategy that preserves this structural element.
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

    # Track molecules with dimethoxybenzene pattern
    molecules_with_dimethoxybenzene = []
    all_nodes = []
    result = False

    def check_dimethoxybenzene(smiles):
        """Check if a molecule contains any dimethoxybenzene pattern"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Check for benzene ring first
            if checker.check_ring_system("benzene", mol):
                if "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("benzene")

            # Count methoxy groups on benzene
            methoxy_count = 0
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic() and atom.GetSymbol() == "C":
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == "O" and neighbor.GetDegree() == 2:
                            for nn in neighbor.GetNeighbors():
                                if (
                                    nn.GetSymbol() == "C"
                                    and nn.GetDegree() == 1
                                    and nn.GetNumImplicitHs() == 3
                                    and nn != atom
                                ):
                                    methoxy_count += 1

            if methoxy_count >= 2:
                if "dimethoxybenzene" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("dimethoxybenzene")
                return True
            return False
        except Exception as e:
            print(f"Error checking dimethoxybenzene pattern: {e}")
            return False

    def collect_all_nodes(node):
        """Collect all nodes in the route"""
        all_nodes.append(node)
        for child in node.get("children", []):
            collect_all_nodes(child)

    def dfs_traverse(node, depth=0):
        """Traverse the route and check for dimethoxybenzene patterns"""
        nonlocal molecules_with_dimethoxybenzene

        if node["type"] == "mol" and "smiles" in node:
            try:
                smiles = node["smiles"]
                if check_dimethoxybenzene(smiles):
                    molecules_with_dimethoxybenzene.append(node)
                    # print(f"Dimethoxybenzene motif detected in: {smiles}") # Original print statement
            except Exception as e:
                print(f"Error processing molecule SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    # First collect all nodes
    collect_all_nodes(route)

    # Then start traversal from the root
    dfs_traverse(route)

    # Check if the final product (root node) has dimethoxybenzene
    final_product_has_dimethoxy = False
    if route["type"] == "mol" and check_dimethoxybenzene(route["smiles"]):
        final_product_has_dimethoxy = True
        # print(f"Final product has dimethoxybenzene motif: {route['smiles']}") # Original print statement
        if {"type": "positional", "details": {"target": "dimethoxybenzene", "position": "final_product"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "dimethoxybenzene", "position": "final_product"}})

    # Check if any starting material has dimethoxybenzene
    starting_materials_with_dimethoxy = []
    starting_material_has_dimethoxy = False
    for node in all_nodes:
        if (
            node["type"] == "mol"
            and node.get("in_stock", False)
            and check_dimethoxybenzene(node["smiles"])
        ):
            starting_materials_with_dimethoxy.append(node["smiles"])
            starting_material_has_dimethoxy = True
            # print(f"Starting material with dimethoxybenzene: {node['smiles']}") # Original print statement

    if starting_material_has_dimethoxy:
        if {"type": "positional", "details": {"target": "dimethoxybenzene", "position": "starting_material"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "dimethoxybenzene", "position": "starting_material"}})

    # Check if dimethoxybenzene is preserved throughout the synthesis
    if final_product_has_dimethoxy and starting_materials_with_dimethoxy:
        # print("Dimethoxybenzene motif is preserved throughout the synthesis") # Original print statement
        result = True

    return result, findings_json
