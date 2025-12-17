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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


ALKYLATION_REACTION_TYPES = [
    "Alkylation of amines",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "S-alkylation of thiols",
    "S-alkylation of thiols (ethyl)",
    "S-alkylation of thiols with alcohols",
    "S-alkylation of thiols with alcohols (ethyl)",
    "Friedel-Crafts alkylation",
    "Friedel-Crafts alkylation with halide",
    "Williamson Ether Synthesis",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that involves the formation of a quaternary carbon center
    through sequential alkylations.
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

    # Track quaternary carbon formation and alkylation reactions
    alkylation_reactions = []
    quaternary_carbon_formed = False
    result = False # Initialize the main boolean result

    def is_quaternary_carbon(atom):
        """Check if a carbon atom is quaternary (connected to 4 other atoms)"""
        if atom.GetAtomicNum() != 6:  # Not carbon
            return False

        # Count non-hydrogen neighbors (explicit and implicit)
        heavy_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)
        return heavy_neighbors == 4

    def get_quaternary_carbons(mol):
        """Get all quaternary carbon atoms in a molecule with their map numbers"""
        result = {}
        for atom in mol.GetAtoms():
            if is_quaternary_carbon(atom):
                map_num = atom.GetAtomMapNum()
                if map_num > 0:  # Only track mapped atoms
                    result[map_num] = atom
        return result

    def is_alkylation_reaction(rsmi):
        """Check if the reaction is an alkylation"""
        for rxn_type in ALKYLATION_REACTION_TYPES:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def dfs_traverse(node, depth=0, quat_carbon_history=None):
        nonlocal quaternary_carbon_formed, alkylation_reactions, findings_json

        if quat_carbon_history is None:
            quat_carbon_history = {}

        if node["type"] == "reaction":
            try:
                # Extract reactants and product from atom-mapped reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                is_current_reaction_alkylation = is_alkylation_reaction(rsmi)

                # Check if this is an alkylation reaction
                if is_current_reaction_alkylation:
                    alkylation_reactions.append((depth, rsmi))
                    print(f"Alkylation reaction detected at depth {depth}: {rsmi}")

                # Convert to RDKit molecules
                reactants_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(products_part)

                if reactants_mol and product_mol:
                    # Get quaternary carbons in reactants and product
                    reactant_quat_carbons = get_quaternary_carbons(reactants_mol)
                    product_quat_carbons = get_quaternary_carbons(product_mol)

                    # Find newly formed quaternary carbons (in product but not in reactants)
                    new_quat_carbon_maps = set(product_quat_carbons.keys()) - set(
                        reactant_quat_carbons.keys()
                    )

                    if new_quat_carbon_maps:
                        # If a new quaternary carbon is formed
                        if "quaternary_carbon_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("quaternary_carbon_formation")

                        if is_current_reaction_alkylation:
                            print(f"Quaternary carbon formation detected at depth {depth}")
                            print(f"New quaternary carbon map numbers: {new_quat_carbon_maps}")

                            # Update quaternary carbon history
                            for map_num in new_quat_carbon_maps:
                                if map_num in quat_carbon_history:
                                    quat_carbon_history[map_num].append(depth)
                                else:
                                    quat_carbon_history[map_num] = [depth]

                            quaternary_carbon_formed = True
                            # Record co-occurrence constraint if both conditions met in this step
                            co_occurrence_constraint = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "alkylation_reaction",
                                        "quaternary_carbon_formation"
                                    ],
                                    "description": "Requires that at least one reaction step in the route is simultaneously an alkylation reaction and results in the formation of a new quaternary carbon."
                                }
                            }
                            if co_occurrence_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(co_occurrence_constraint)

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # For molecule nodes, pass the quaternary carbon history to children
        elif node["type"] == "mol" and not node.get("in_stock", False):
            # Only process non-starting materials
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                quat_carbons = get_quaternary_carbons(mol)
                # Update history with current quaternary carbons
                for map_num, atom in quat_carbons.items():
                    if map_num not in quat_carbon_history:
                        quat_carbon_history[map_num] = []

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth, quat_carbon_history.copy())

    # Start traversal
    dfs_traverse(route)

    # Check if we have sequential alkylations (at least 2)
    sequential_alkylations = len(alkylation_reactions) >= 2

    if sequential_alkylations:
        count_constraint = {
            "type": "count",
            "details": {
                "target": "alkylation_reaction",
                "operator": ">=",
                "value": 2,
                "description": "Checks for at least two alkylation reactions anywhere in the synthesis route. 'alkylation_reaction' refers to any of the specific alkylation types listed in atomic_checks."
            }
        }
        if count_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(count_constraint)

    if quaternary_carbon_formed and sequential_alkylations:
        print(f"Found quaternary carbon formation through sequential alkylations")
        print(f"Alkylation reactions: {alkylation_reactions}")
        result = True

    return result, findings_json
