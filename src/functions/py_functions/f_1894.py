#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a convergent synthesis strategy involving amide coupling
    with a trifluoromethyl-containing aromatic fragment.
    """
    # Track if we found the required elements
    amide_coupling_found = False
    trifluoromethyl_fragment_involved = False

    # Track reaction nodes to analyze convergence
    reaction_nodes = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal amide_coupling_found, trifluoromethyl_fragment_involved

        if path is None:
            path = []

        current_path = path + [node]

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reaction_nodes.append((node, depth, current_path))

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide coupling reaction
            is_amide_coupling = False
            for reaction_type in [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Schotten-Baumann_amide",
            ]:
                if checker.check_reaction(reaction_type, rsmi):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling reaction: {reaction_type} at depth {depth}"
                    )
                    break

            # If not detected by reaction type, check for functional group transformation
            if not is_amide_coupling:
                # Check if reactants contain carboxylic acid and amine
                has_carboxylic_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants if r
                )
                has_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants if r
                )

                # Check if product contains amide
                has_amide = checker.check_fg(
                    "Primary amide", product
                ) or checker.check_fg("Secondary amide", product)

                if has_carboxylic_acid and has_amine and has_amide:
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling by functional group analysis at depth {depth}"
                    )

            if is_amide_coupling:
                amide_coupling_found = True

                # Check if any reactant contains trifluoromethyl group on aromatic ring
                for reactant in reactants:
                    if reactant and checker.check_fg("Trifluoro group", reactant):
                        # Verify it's on an aromatic fragment
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Check if the molecule has aromatic atoms
                            has_aromatic = False
                            for atom in mol.GetAtoms():
                                if atom.GetIsAromatic():
                                    has_aromatic = True
                                    break

                            if has_aromatic:
                                trifluoromethyl_fragment_involved = True
                                print(
                                    f"Trifluoromethyl aromatic fragment involved in amide coupling at depth {depth}"
                                )
                                break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Traverse the route to collect all reaction nodes
    dfs_traverse(route)

    # Check if the synthesis is convergent (multiple branches coming together)
    is_convergent = False
    for node, depth, path in reaction_nodes:
        if amide_coupling_found and trifluoromethyl_fragment_involved:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # A convergent synthesis typically has multiple reactants
            if len(reactants) >= 2:
                is_convergent = True
                print(
                    f"Convergent synthesis detected with {len(reactants)} reactants at depth {depth}"
                )
                break

    result = (
        amide_coupling_found and trifluoromethyl_fragment_involved and is_convergent
    )
    print(f"Final result: {result}")
    print(f"- Amide coupling found: {amide_coupling_found}")
    print(
        f"- Trifluoromethyl aromatic fragment involved: {trifluoromethyl_fragment_involved}"
    )
    print(f"- Convergent synthesis: {is_convergent}")

    return result
