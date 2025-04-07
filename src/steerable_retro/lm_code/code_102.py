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
    Detects if the synthesis route involves C-N bond formation in the late stage (depth 0-1).
    """
    late_stage_cn_formation = False

    # List of reaction types that typically form C-N bonds
    cn_bond_forming_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Aminolysis of esters",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Ugi reaction",
        "Schotten-Baumann to ester",
        "Schotten-Baumann_amide",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_formation

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a C-N bond forming reaction using the checker
                for reaction_type in cn_bond_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Late-stage C-N bond formation detected: {reaction_type} at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        late_stage_cn_formation = True
                        return

                # If no specific reaction type matched, check for C-N bond formation directly
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(reactant_mols):
                    # Check for new C-N bonds using atom mapping
                    # Get all C-N bonds in the product
                    product_cn_bonds = []
                    for bond in product_mol.GetBonds():
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 7) or (
                            begin_atom.GetAtomicNum() == 7 and end_atom.GetAtomicNum() == 6
                        ):
                            # Get atom maps if available
                            begin_map = (
                                begin_atom.GetProp("molAtomMapNumber")
                                if begin_atom.HasProp("molAtomMapNumber")
                                else None
                            )
                            end_map = (
                                end_atom.GetProp("molAtomMapNumber")
                                if end_atom.HasProp("molAtomMapNumber")
                                else None
                            )
                            if begin_map and end_map:
                                product_cn_bonds.append((begin_map, end_map))

                    # Check if these C-N bonds exist in reactants
                    for begin_map, end_map in product_cn_bonds:
                        bond_exists_in_reactants = False
                        for r_mol in reactant_mols:
                            # Find atoms with these mapping numbers
                            begin_atom_idx = None
                            end_atom_idx = None
                            for atom in r_mol.GetAtoms():
                                if atom.HasProp("molAtomMapNumber"):
                                    atom_map = atom.GetProp("molAtomMapNumber")
                                    if atom_map == begin_map:
                                        begin_atom_idx = atom.GetIdx()
                                    elif atom_map == end_map:
                                        end_atom_idx = atom.GetIdx()

                            # Check if there's a bond between these atoms
                            if begin_atom_idx is not None and end_atom_idx is not None:
                                bond = r_mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
                                if bond:
                                    bond_exists_in_reactants = True
                                    break

                        if not bond_exists_in_reactants:
                            print(
                                f"Late-stage C-N bond formation detected through atom mapping at depth {depth}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            late_stage_cn_formation = True
                            return

                # Check for specific functional groups that involve C-N bonds
                if not late_stage_cn_formation:
                    # Check if product has amines, amides, etc. that reactants don't have
                    cn_functional_groups = [
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                        "Aniline",
                        "Urea",
                        "Carbamate",
                    ]

                    for fg in cn_functional_groups:
                        if checker.check_fg(fg, product_smiles):
                            # Check if this FG exists in any reactant
                            fg_exists_in_reactants = any(
                                checker.check_fg(fg, r) for r in reactants_smiles
                            )
                            if not fg_exists_in_reactants:
                                print(
                                    f"Late-stage C-N bond formation detected: new {fg} at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                late_stage_cn_formation = True
                                return

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            if (
                not late_stage_cn_formation
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_cn_formation
