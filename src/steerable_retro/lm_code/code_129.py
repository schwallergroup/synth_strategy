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
    Detects a linear synthetic strategy where a single aromatic position
    is sequentially functionalized without branching or convergent steps.
    """
    # Track the number of sequential aromatic functionalizations
    aromatic_functionalizations = 0
    max_depth = 0

    # List to store the reactions in order (for checking sequential nature)
    reactions_list = []

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_functionalizations, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Try to create RDKit molecules
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Skip if any molecule failed to parse
                if not product_mol or None in reactant_mols:
                    return

                # Check for aromatic functionalization reactions
                is_functionalization = False

                # Check for common aromatic functionalization reactions
                if (
                    checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                    or checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Friedel-Crafts alkylation", rsmi)
                    or checker.check_reaction("Friedel-Crafts acylation", rsmi)
                    or checker.check_reaction("Aromatic hydroxylation", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Acylation of olefines by aldehydes", rsmi)
                ):

                    # Find the main reactant (the one with an aromatic ring that's being modified)
                    main_reactant = None
                    main_reactant_idx = -1

                    for i, reactant in enumerate(reactant_mols):
                        if reactant.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")):
                            # Check if this reactant shares a significant substructure with the product
                            mcs = rdFMCS.FindMCS(
                                [reactant, product_mol],
                                completeRingsOnly=True,
                                ringMatchesRingOnly=True,
                            )
                            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

                            if (
                                mcs_mol
                                and reactant.HasSubstructMatch(mcs_mol)
                                and product_mol.HasSubstructMatch(mcs_mol)
                            ):
                                # If the MCS contains an aromatic ring, this is likely our main reactant
                                if mcs_mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")):
                                    main_reactant = reactant
                                    main_reactant_idx = i
                                    break

                    if main_reactant and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("c1ccccc1")
                    ):
                        is_functionalization = True
                        print(f"Found aromatic functionalization reaction at depth {depth}: {rsmi}")

                # Also check for functional group changes on the aromatic ring
                elif any(
                    reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1"))
                    for reactant_mol in reactant_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")):
                    # Check for specific functional group changes
                    fg_changes = [
                        ("Aromatic halide", "Aromatic halide"),
                        ("Nitro group", "Nitro group"),
                        ("Phenol", "Phenol"),
                        ("Aniline", "Aniline"),
                        ("Aromatic thiol", "Aromatic thiol"),
                        ("Nitrile", "Nitrile"),
                        ("Carboxylic acid", "Carboxylic acid"),
                        ("Ester", "Ester"),
                        ("Aldehyde", "Aldehyde"),
                        ("Ketone", "Ketone"),
                    ]

                    for fg_name, _ in fg_changes:
                        # Check if the product has a functional group that at least one reactant doesn't
                        if checker.check_fg(fg_name, product_smiles):
                            has_fg_in_reactants = False
                            for r_smiles in reactants_smiles:
                                if checker.check_fg(fg_name, r_smiles):
                                    has_fg_in_reactants = True
                                    break

                            if not has_fg_in_reactants:
                                is_functionalization = True
                                print(
                                    f"Found aromatic functionalization (FG change: {fg_name}) at depth {depth}"
                                )
                                break

                if is_functionalization:
                    reactions_list.append((depth, rsmi))
                    aromatic_functionalizations += 1

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total aromatic functionalizations found: {aromatic_functionalizations}")
    print(f"Maximum depth in synthesis: {max_depth}")

    # Sort reactions by depth
    reactions_list.sort(key=lambda x: x[0])

    # Print the reactions for debugging
    for depth, rsmi in reactions_list:
        print(f"Depth {depth}: {rsmi}")

    # Check if we have enough functionalizations
    if aromatic_functionalizations < 2:
        print("Not enough aromatic functionalizations (need at least 2)")
        return False

    # Check if the reactions form a reasonably continuous sequence
    # Allow for some gaps in the sequence (e.g., protection/deprotection steps)
    is_continuous = True
    if len(reactions_list) >= 2:
        depths = [r[0] for r in reactions_list]
        for i in range(1, len(depths)):
            # Allow gaps of up to 4 in depth (to account for intermediate steps)
            if depths[i] - depths[i - 1] > 4:
                print(f"Gap in sequence between depths {depths[i-1]} and {depths[i]}")
                is_continuous = False
                break
    else:
        is_continuous = False

    # Check if the synthesis is predominantly linear
    is_linear = aromatic_functionalizations >= max(1, max_depth // 3)

    print(f"Is continuous: {is_continuous}")
    print(f"Is linear: {is_linear}")

    # Strategy is present if:
    # 1. We have at least 2 aromatic functionalizations
    # 2. The functionalizations form a reasonably continuous sequence
    # 3. The synthesis has a significant linear component
    return aromatic_functionalizations >= 2 and is_continuous and is_linear
