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
    This function detects convergent synthesis where multiple fragments
    are combined in late-stage reactions.
    """
    fragment_count = 0
    convergent = False
    unique_fragments = set()  # Track unique fragments
    coupling_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_count, convergent, coupling_reactions

        if node["type"] == "reaction" and depth <= 3:  # Focus on late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                is_coupling = False
                coupling_rxn_types = [
                    "Suzuki",
                    "Negishi",
                    "Stille",
                    "Heck",
                    "Sonogashira",
                    "Buchwald-Hartwig",
                    "Ullmann",
                    "Kumada",
                ]

                for rxn_type in coupling_rxn_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        coupling_reactions += 1
                        print(f"Found {rxn_type} coupling reaction at depth {depth}")
                        break

                # Count distinct reactants with significant complexity
                significant_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Consider size, rings, or functional group complexity
                            if (
                                mol.GetNumAtoms() > 6
                                or mol.GetNumRings() > 0
                                or any(
                                    checker.check_fg(fg, reactant)
                                    for fg in [
                                        "Aromatic halide",
                                        "Boronic acid",
                                        "Boronic ester",
                                        "Alkyne",
                                    ]
                                )
                            ):

                                # Add canonical SMILES to avoid double counting
                                canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                                if canonical_smiles not in unique_fragments:
                                    unique_fragments.add(canonical_smiles)
                                    significant_reactants += 1
                                    print(
                                        f"Found significant fragment: {canonical_smiles} at depth {depth}"
                                    )
                    except Exception as e:
                        print(f"Error processing reactant: {e}")
                        continue

                if significant_reactants >= 2:
                    fragment_count += significant_reactants
                    print(f"Found {significant_reactants} significant fragments at depth {depth}")

                    # Check if product is significantly more complex than individual reactants
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.GetNumAtoms() > max(
                            [
                                Chem.MolFromSmiles(r).GetNumAtoms()
                                for r in reactants
                                if Chem.MolFromSmiles(r)
                            ]
                        ):
                            print(
                                f"Product is more complex than individual reactants at depth {depth}"
                            )
                            convergent = True
                    except Exception as e:
                        print(f"Error comparing product complexity: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Consider convergent if multiple significant fragments are combined
    # or if coupling reactions are used in late stages
    if fragment_count >= 3 or (fragment_count >= 2 and coupling_reactions >= 1):
        convergent = True
        print(
            f"Convergent synthesis detected with {fragment_count} fragments and {coupling_reactions} coupling reactions"
        )

    return convergent
