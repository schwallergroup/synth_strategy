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
    Detects if the synthetic route forms an oxadiazolone ring via cyclization of a hydrazide.
    """
    hydrazide_info = {"detected": False, "smiles": None, "depth": -1}
    oxadiazolone_info = {"detected": False, "smiles": None, "depth": -1}
    direct_conversion = False

    def is_oxadiazolone(mol_smiles):
        """Check if the molecule contains an oxadiazolone (1,3,4-oxadiazol-5(2H)-one) structure"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return False

        # Various patterns for 1,3,4-oxadiazol-5(2H)-one in different forms
        patterns = [
            Chem.MolFromSmarts("c1nnoc1=O"),
            Chem.MolFromSmarts("O=c1[nH]nco1"),
            Chem.MolFromSmarts("O=c1onc[nH]1"),
            Chem.MolFromSmarts("O=C1ONC=N1"),
            Chem.MolFromSmarts("O=c1nc(O)on1"),
            Chem.MolFromSmarts("O=C1N=CON1"),
            Chem.MolFromSmarts("O=c1onc(O)n1"),
        ]

        for pattern in patterns:
            if pattern and mol.HasSubstructMatch(pattern):
                print(f"Oxadiazolone pattern match found in: {mol_smiles}")
                return True

        # Check if it has both oxadiazole ring and carbonyl group
        if checker.check_ring("oxadiazole", mol_smiles):
            if checker.check_fg("Ketone", mol_smiles) or checker.check_fg("Amide", mol_smiles):
                print(f"Oxadiazole with carbonyl group found in: {mol_smiles}")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal hydrazide_info, oxadiazolone_info, direct_conversion

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for hydrazide in molecule nodes
            if checker.check_fg("Acylhydrazine", mol_smiles) and (
                not hydrazide_info["detected"] or depth > hydrazide_info["depth"]
            ):
                hydrazide_info["detected"] = True
                hydrazide_info["smiles"] = mol_smiles
                hydrazide_info["depth"] = depth
                print(f"Hydrazide detected at depth {depth}: {mol_smiles}")

            # Check for oxadiazolone in molecule nodes
            if is_oxadiazolone(mol_smiles) and (
                not oxadiazolone_info["detected"] or depth < oxadiazolone_info["depth"]
            ):
                oxadiazolone_info["detected"] = True
                oxadiazolone_info["smiles"] = mol_smiles
                oxadiazolone_info["depth"] = depth
                print(f"Oxadiazolone ring detected at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydrazide formation reaction
                if checker.check_fg("Acylhydrazine", product) and (
                    not hydrazide_info["detected"] or depth > hydrazide_info["depth"]
                ):
                    for reactant in reactants:
                        if not checker.check_fg("Acylhydrazine", reactant):
                            hydrazide_info["detected"] = True
                            hydrazide_info["smiles"] = product
                            hydrazide_info["depth"] = depth
                            print(f"Hydrazide formation detected at depth {depth}")
                            break

                # Check for oxadiazolone formation reaction
                if is_oxadiazolone(product) and (
                    not oxadiazolone_info["detected"] or depth < oxadiazolone_info["depth"]
                ):
                    oxadiazolone_info["detected"] = True
                    oxadiazolone_info["smiles"] = product
                    oxadiazolone_info["depth"] = depth
                    print(f"Oxadiazolone ring detected in product at depth {depth}: {product}")

                    # Check if any reactant contains hydrazide for direct conversion
                    for reactant in reactants:
                        if checker.check_fg("Acylhydrazine", reactant):
                            direct_conversion = True
                            print(
                                f"Direct conversion: Hydrazide to oxadiazolone detected at depth {depth}"
                            )
                            break

                    # Additional reaction checks for cyclization
                    if not direct_conversion:
                        reactants_combined = ".".join(reactants)
                        if checker.check_fg("Acylhydrazine", reactants_combined):
                            # Check for cyclization/dehydration reactions
                            cyclization_reactions = [
                                "Formation of NOS Heterocycles",
                                "1,2,4-oxadiazol-5(2H)-one synthesis from nitrile, hydrogen carbonate, and hydroxylamine",
                                "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                                "{oxadiazole}",
                            ]

                            for rxn_type in cyclization_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    direct_conversion = True
                                    print(
                                        f"Cyclization reaction {rxn_type} detected at depth {depth}"
                                    )
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check for correct sequence in retrosynthetic traversal
    # In retrosynthesis, lower depth means later stage (product)
    correct_sequence = False
    if hydrazide_info["detected"] and oxadiazolone_info["detected"]:
        if oxadiazolone_info["depth"] < hydrazide_info["depth"]:
            correct_sequence = True
            print(
                f"Correct sequence detected: oxadiazolone at depth {oxadiazolone_info['depth']} <- hydrazide at depth {hydrazide_info['depth']}"
            )

    # For debugging
    print(f"Hydrazide detected: {hydrazide_info['detected']} at depth {hydrazide_info['depth']}")
    print(
        f"Oxadiazolone detected: {oxadiazolone_info['detected']} at depth {oxadiazolone_info['depth']}"
    )
    print(f"Direct conversion: {direct_conversion}")
    print(f"Correct sequence: {correct_sequence}")

    # Return true if we found both a hydrazide and an oxadiazolone in the correct sequence
    # OR if we found a direct conversion from hydrazide to oxadiazolone
    return hydrazide_info["detected"] and (
        (oxadiazolone_info["detected"] and correct_sequence) or direct_conversion
    )
