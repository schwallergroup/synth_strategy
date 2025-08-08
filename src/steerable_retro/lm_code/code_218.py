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
    This function detects a synthetic strategy involving isocyanate intermediates
    that undergo cyclization to form nitrogen heterocycles.
    """
    # Track molecules with isocyanates and their positions in the route
    isocyanate_molecules = set()
    # Track cyclization reactions that form nitrogen heterocycles
    cyclization_reactions = []
    # Track the depth of each node to determine sequence
    node_depths = {}

    def dfs_traverse(node, depth=0):
        # Store the depth of this node
        if node["type"] == "mol" and "smiles" in node and node["smiles"]:
            node_depths[node["smiles"]] = depth

            # Check for isocyanate in molecules
            if checker.check_fg("Isocyanate", node["smiles"]):
                print(f"Found isocyanate in molecule: {node['smiles']} at depth {depth}")
                isocyanate_molecules.add(node["smiles"])

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a cyclization reaction forming a nitrogen heterocycle
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Count rings in reactants and product
                reactant_ring_count = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)
                product_ring_count = product_mol.GetRingInfo().NumRings()

                # Check if product has more rings than reactants (cyclization)
                if product_ring_count > reactant_ring_count:
                    # Check if the new ring contains nitrogen (nitrogen heterocycle)
                    nitrogen_in_ring = False
                    ring_info = product_mol.GetRingInfo()
                    for ring_atoms in ring_info.AtomRings():
                        for atom_idx in ring_atoms:
                            atom = product_mol.GetAtomWithIdx(atom_idx)
                            if atom.GetAtomicNum() == 7:  # Nitrogen
                                nitrogen_in_ring = True
                                break
                        if nitrogen_in_ring:
                            break

                    if nitrogen_in_ring:
                        print(f"Found cyclization forming nitrogen heterocycle: {rsmi}")
                        # Store the reaction and its reactants for later analysis
                        cyclization_reactions.append((rsmi, reactants, product, depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both isocyanates and cyclization reactions
    if not isocyanate_molecules or not cyclization_reactions:
        return False

    # Check if any cyclization reaction uses an isocyanate
    for rxn, reactants, product, rxn_depth in cyclization_reactions:
        for reactant in reactants:
            # If this reactant is or contains an isocyanate
            for isocyanate_mol in isocyanate_molecules:
                isocyanate_depth = node_depths.get(isocyanate_mol, -1)
                # Check if the isocyanate appears before the cyclization in the synthesis route
                # (higher depth means earlier in the synthesis)
                if isocyanate_depth >= rxn_depth and (
                    reactant == isocyanate_mol or checker.check_fg("Isocyanate", reactant)
                ):
                    print(
                        f"Confirmed isocyanate cyclization strategy: isocyanate at depth {isocyanate_depth} used in cyclization at depth {rxn_depth}"
                    )
                    return True

    return False
