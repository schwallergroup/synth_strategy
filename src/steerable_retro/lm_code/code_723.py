#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
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

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects preservation of cyclopropyl motifs throughout the synthesis.
    Returns True if the final product contains a cyclopropyl ring AND all intermediates
    in the main synthetic pathway leading to it also contain the cyclopropyl ring.
    """
    has_cyclopropyl_in_product = False
    cyclopropyl_atoms_in_product = set()

    def get_atom_mapping(reaction_smiles):
        """Extract atom mapping from reaction SMILES"""
        try:
            rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            product = rxn.GetProducts()[0]
            reactants = rxn.GetReactants()

            # Get atom mapping dictionary (product atom idx -> reactant idx, reactant atom idx)
            atom_mapping = {}
            for atom in product.GetAtoms():
                if atom.HasProp("molAtomMapNumber"):
                    map_num = atom.GetProp("molAtomMapNumber")
                    for r_idx, reactant in enumerate(reactants):
                        for a in reactant.GetAtoms():
                            if (
                                a.HasProp("molAtomMapNumber")
                                and a.GetProp("molAtomMapNumber") == map_num
                            ):
                                atom_mapping[atom.GetIdx()] = (r_idx, a.GetIdx())
            return atom_mapping
        except Exception as e:
            print(f"Error in get_atom_mapping: {e}")
            return {}

    def get_cyclopropyl_atoms(mol_smiles):
        """Get atom indices involved in cyclopropyl rings"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return set()

            # Use RingInfo to find 3-membered carbon rings
            ring_info = mol.GetRingInfo()
            cyclopropyl_atoms = set()

            for ring in ring_info.AtomRings():
                if len(ring) == 3:
                    # Check if all atoms in the ring are carbon
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                        cyclopropyl_atoms.update(ring)

            return cyclopropyl_atoms
        except Exception as e:
            print(f"Error in get_cyclopropyl_atoms: {e}")
            return set()

    def dfs_traverse(node, depth=0, tracked_cyclopropyl=None):
        nonlocal has_cyclopropyl_in_product, cyclopropyl_atoms_in_product

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cyclopropyl = checker.check_ring("cyclopropane", mol_smiles)

            if depth == 0:  # Final product
                has_cyclopropyl_in_product = has_cyclopropyl
                if has_cyclopropyl:
                    cyclopropyl_atoms_in_product = get_cyclopropyl_atoms(mol_smiles)
                    print(
                        f"Final product contains cyclopropyl ring with atoms: {cyclopropyl_atoms_in_product}"
                    )
                    return True
                else:
                    print("Final product does not contain cyclopropyl ring")
                    return False

            # For intermediates, only check if we're tracking cyclopropyl
            if tracked_cyclopropyl is not None:
                if not has_cyclopropyl:
                    print(f"Cyclopropyl not preserved at depth {depth}, molecule: {mol_smiles}")
                    return False
                print(f"Cyclopropyl preserved at depth {depth}, molecule: {mol_smiles}")

            return True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # For reaction nodes, determine which reactants should have cyclopropyl
            # based on atom mapping from the product
            if tracked_cyclopropyl is not None:
                atom_mapping = get_atom_mapping(rsmi)

                # Determine which reactants contain atoms from the cyclopropyl ring
                reactants_with_cyclopropyl = set()
                for atom_idx in tracked_cyclopropyl:
                    if atom_idx in atom_mapping:
                        reactants_with_cyclopropyl.add(atom_mapping[atom_idx][0])

                print(f"Reactants with cyclopropyl at depth {depth}: {reactants_with_cyclopropyl}")

                # Update tracked_cyclopropyl for each child based on atom mapping
                for i, child in enumerate(node.get("children", [])):
                    if child["type"] == "mol":
                        if i in reactants_with_cyclopropyl:
                            # This reactant should have cyclopropyl
                            if not dfs_traverse(child, depth + 1, True):
                                return False
                        else:
                            # This reactant doesn't need to have cyclopropyl
                            if not dfs_traverse(child, depth + 1, None):
                                return False
                return True

        # Continue traversal for all children
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1, tracked_cyclopropyl):
                return False

        return True

    # Start traversal from the root
    result = dfs_traverse(route)

    # Return True only if product has cyclopropyl AND it's preserved throughout
    return has_cyclopropyl_in_product and result
