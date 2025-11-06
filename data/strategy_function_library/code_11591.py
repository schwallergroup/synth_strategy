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
    This function detects a synthetic strategy involving the formation of a quaternary carbon center
    through alpha-functionalization of a cyclic ketone.
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

    quaternary_center_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal quaternary_center_formed, findings_json

        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if any reactant is a cyclic ketone
                for reactant_smiles in reactants_smiles:
                    if not checker.check_fg("Ketone", reactant_smiles):
                        continue
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if not reactant_mol:
                        continue

                    # Get ring information
                    ring_info = reactant_mol.GetRingInfo()
                    if ring_info.NumRings() == 0:
                        continue  # Skip if no rings

                    # Get ketone functional group atom indices
                    ketone_atoms = checker.get_fg_atom_indices("Ketone", reactant_smiles)
                    if not ketone_atoms:
                        continue

                    # Parse the product molecule
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if not product_mol:
                        continue

                    # Check each ketone group
                    for atom_indices_tuple in ketone_atoms:
                        # Flatten the tuple of tuples if necessary
                        atom_indices = (
                            atom_indices_tuple[0]
                            if isinstance(atom_indices_tuple[0], tuple)
                            else atom_indices_tuple
                        )

                        # Find the carbonyl carbon (connected to oxygen)
                        carbonyl_carbon_idx = None
                        for idx in atom_indices:
                            atom = reactant_mol.GetAtomWithIdx(idx)
                            if atom.GetSymbol() == "C":
                                for neighbor in atom.GetNeighbors():
                                    if (
                                        neighbor.GetSymbol() == "O"
                                        and neighbor.GetIdx() in atom_indices
                                    ):
                                        carbonyl_carbon_idx = idx
                                        break
                                if carbonyl_carbon_idx is not None:
                                    break

                        if carbonyl_carbon_idx is None:
                            continue

                        # Find alpha carbons (adjacent to carbonyl carbon and in the same ring)
                        alpha_carbons = []
                        carbonyl_atom = reactant_mol.GetAtomWithIdx(carbonyl_carbon_idx)
                        for neighbor in carbonyl_atom.GetNeighbors():
                            if (
                                neighbor.GetSymbol() == "C"
                                and ring_info.AreAtomsInSameRing(
                                    carbonyl_carbon_idx, neighbor.GetIdx()
                                )
                            ):
                                alpha_carbons.append(neighbor.GetIdx())

                        if not alpha_carbons:
                            continue

                        # Check if any alpha carbon is tertiary in reactant and becomes quaternary in product
                        for alpha_idx in alpha_carbons:
                            alpha_atom = reactant_mol.GetAtomWithIdx(alpha_idx)

                            # Check if alpha carbon is tertiary (connected to 3 atoms)
                            if alpha_atom.GetDegree() == 3:
                                # Get atom mapping number to track this atom in the product
                                alpha_map_num = (
                                    alpha_atom.GetProp("molAtomMapNumber")
                                    if alpha_atom.HasProp("molAtomMapNumber")
                                    else None
                                )

                                if alpha_map_num:
                                    # Find the corresponding atom in the product
                                    for product_atom in product_mol.GetAtoms():
                                        if (
                                            product_atom.HasProp("molAtomMapNumber")
                                            and product_atom.GetProp("molAtomMapNumber")
                                            == alpha_map_num
                                        ):

                                            # Check if it's quaternary in the product (connected to 4 atoms)
                                            if product_atom.GetDegree() >= 4:
                                                print(
                                                    f"Detected quaternary center formation at depth {depth}"
                                                )
                                                print(f"Reactant: {reactant_smiles}")
                                                print(f"Product: {product_smiles}")
                                                print(
                                                    f"Alpha carbon with map number {alpha_map_num} became quaternary"
                                                )
                                                quaternary_center_formed = True
                                                findings_json["atomic_checks"]["named_reactions"].append("quaternary_center_formation")
                                                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "quaternary_center_formation", "operator": ">=", "value": 1}})
                                                return  # Exit once found

            except Exception as e:
                print(f"Error in quaternary center detection: {str(e)}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node.get("type") != "reaction":  # If current node is 'chemical'
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)
            if quaternary_center_formed:
                return  # Exit early if pattern found

    # Start traversal from the root
    dfs_traverse(route)

    return quaternary_center_formed, findings_json
