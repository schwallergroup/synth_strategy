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
    Detects a strategy where a nitrile group is transformed as part of heterocyclic ring formation
    """
    found_nitrile_transformation = False

    # List of nitrogen heterocycles to check for
    n_heterocycles = [
        "tetrazole",
        "triazole",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of nitrile-to-heterocycle reactions
    nitrile_heterocycle_reactions = [
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "A3 coupling to imidazoles",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_transformation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for specific nitrile-to-heterocycle reactions
                for reaction_type in nitrile_heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found {reaction_type} at depth {depth}")
                        found_nitrile_transformation = True
                        return

                # Check for nitrile in reactants
                nitrile_reactant_idx = None
                for idx, reactant_smiles in enumerate(reactants_smiles):
                    if checker.check_fg("Nitrile", reactant_smiles):
                        nitrile_reactant_idx = idx
                        break

                if nitrile_reactant_idx is not None:
                    # Convert to RDKit molecules
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                    if product_mol is None or None in reactant_mols:
                        return

                    # Count rings in product and reactants
                    product_ring_count = product_mol.GetRingInfo().NumRings()
                    reactant_ring_counts = [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    total_reactant_rings = sum(reactant_ring_counts)

                    # Check if product has more rings (ring formation occurred)
                    if product_ring_count > total_reactant_rings:
                        # Check if the product contains a nitrogen heterocycle
                        for heterocycle in n_heterocycles:
                            if checker.check_ring(heterocycle, product_smiles):
                                print(
                                    f"Found nitrile transformation to {heterocycle} at depth {depth}"
                                )
                                found_nitrile_transformation = True
                                return

                        # Fallback to manual check for nitrogen in rings
                        nitrogen_in_ring = False
                        for ring in product_mol.GetRingInfo().AtomRings():
                            for atom_idx in ring:
                                atom = product_mol.GetAtomWithIdx(atom_idx)
                                if atom.GetAtomicNum() == 7:  # Nitrogen
                                    nitrogen_in_ring = True
                                    break
                            if nitrogen_in_ring:
                                break

                        if nitrogen_in_ring:
                            print(
                                f"Found nitrile transformation in heterocycle formation at depth {depth}"
                            )
                            found_nitrile_transformation = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitrile_transformation
