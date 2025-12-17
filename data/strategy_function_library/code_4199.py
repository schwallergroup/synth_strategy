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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving early formation of a quaternary carbon
    center with a nitrile group.
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

    # Initialize tracking variable
    has_quaternary_nitrile_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_quaternary_nitrile_formation, findings_json

        if node["type"] == "reaction" and depth >= 3:  # Early stage = higher depth (at least 3)
            # Extract reaction information
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if product has nitrile group
                if checker.check_fg("Nitrile", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    # Get the product molecule
                    product_mol = Chem.MolFromSmiles(product)
                    if not product_mol:
                        return

                    # Find quaternary carbons in product
                    quaternary_carbons = []
                    for atom in product_mol.GetAtoms():
                        if atom.GetAtomicNum() == 6:  # Carbon atom
                            # Count non-hydrogen connections
                            non_h_neighbors = [
                                n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1
                            ]
                            if len(non_h_neighbors) == 4:
                                quaternary_carbons.append(atom.GetIdx())

                    if quaternary_carbons:
                        # Check if any quaternary carbon is connected to a nitrile group
                        nitrile_connected = False
                        for qc_idx in quaternary_carbons:
                            qc_atom = product_mol.GetAtomWithIdx(qc_idx)
                            for neighbor in qc_atom.GetNeighbors():
                                # Check if neighbor is carbon and part of a nitrile
                                if neighbor.GetAtomicNum() == 6:
                                    for nn in neighbor.GetNeighbors():
                                        if (
                                            nn.GetAtomicNum() == 7 and nn.GetDegree() == 1
                                        ):  # Nitrogen with one connection
                                            nitrile_connected = True
                                            break

                        if nitrile_connected:
                            # Check if this quaternary carbon was formed in this reaction
                            # by verifying it doesn't exist in reactants
                            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                            reactant_mols = [m for m in reactant_mols if m is not None]

                            # Check if any reactant has a quaternary carbon connected to nitrile
                            reactant_has_quat_nitrile = False
                            for r_mol in reactant_mols:
                                if r_mol and checker.check_fg("Nitrile", Chem.MolToSmiles(r_mol)):
                                    for atom in r_mol.GetAtoms():
                                        if atom.GetAtomicNum() == 6:  # Carbon atom
                                            non_h_neighbors = [
                                                n
                                                for n in atom.GetNeighbors()
                                                if n.GetAtomicNum() != 1
                                            ]
                                            if len(non_h_neighbors) == 4:
                                                # Check if connected to nitrile
                                                for neighbor in atom.GetNeighbors():
                                                    if neighbor.GetAtomicNum() == 6:
                                                        for nn in neighbor.GetNeighbors():
                                                            if (
                                                                nn.GetAtomicNum() == 7
                                                                and nn.GetDegree() == 1
                                                            ):
                                                                reactant_has_quat_nitrile = True
                                                                break

                            if not reactant_has_quat_nitrile:
                                has_quaternary_nitrile_formation = True
                                # Add to findings_json when the condition is met
                                if "quaternary_nitrile_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("quaternary_nitrile_formation")
                                
                                # Add structural constraints if not already present
                                positional_constraint = {"type": "positional", "details": {"target": "quaternary_nitrile_formation", "position": "early_stage"}}
                                if positional_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(positional_constraint)
                                
                                count_constraint = {"type": "count", "details": {"target": "quaternary_nitrile_formation", "operator": ">=", "value": 1}}
                                if count_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(count_constraint)

            except Exception:
                pass

        # Traverse children (depth increases as we go deeper in retrosynthesis)
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for chemical children
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for reaction children
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_quaternary_nitrile_formation, findings_json
