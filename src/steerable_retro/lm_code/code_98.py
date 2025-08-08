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
    Detects if a key C-C bond formation occurs in the final step of the synthesis.
    """
    has_late_stage_cc_bond = False

    # Expanded list of C-C bond forming reactions to check
    cc_bond_forming_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Heck terminal vinyl",
        "Oxidative Heck reaction",
        "Oxidative Heck reaction with vinyl ester",
        "Heck reaction with vinyl ester and amine",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_aryl OTf",
        "Stille reaction_vinyl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Grignard from aldehyde to alcohol",
        "Grignard from ketone to alcohol",
        "Grignard from nitrile to ketone",
        "Wittig reaction with triphenylphosphorane",
        "Wittig with Phosphonium",
        "Aldol condensation",
        "Michael addition",
        "Michael addition methyl",
        "Diels-Alder",
        "Friedel-Crafts alkylation",
        "Friedel-Crafts alkylation with halide",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl OTf",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_acyl halide",
        "Sonogashira acetylene_acyl halide",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Knoevenagel Condensation",
        "Catellani reaction ortho",
        "Catellani reaction para",
        "beta C(sp3) arylation",
        "decarboxylative_coupling",
        "A3 coupling",
        "A3 coupling to imidazoles",
        "Alkyne-imine cycloaddition",
        "Pauson-Khand reaction",
        "Julia Olefination",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cc_bond

        # Check if this is a late-stage reaction (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a known C-C bond forming reaction
            for rxn_type in cc_bond_forming_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected C-C bond forming reaction: {rxn_type}")
                    has_late_stage_cc_bond = True
                    return

            # If not a known reaction type, check by analyzing atom-mapped reaction
            try:
                # Split reaction into reactants and product
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                reactant_mols = [m for m in reactant_mols if m is not None]
                product_mol = Chem.MolFromSmiles(product_part)

                if not product_mol or not reactant_mols:
                    print("Failed to create molecule objects")
                    return

                # Get atom maps from product
                product_c_atoms = {}
                for atom in product_mol.GetAtoms():
                    if (
                        atom.GetAtomicNum() == 6 and atom.GetAtomMapNum() > 0
                    ):  # Carbon atoms with map
                        product_c_atoms[atom.GetAtomMapNum()] = atom.GetIdx()

                # Check for new C-C bonds in product that weren't in reactants
                new_cc_bond_found = False

                # For each pair of mapped carbon atoms in product
                for map_num1, idx1 in product_c_atoms.items():
                    for map_num2, idx2 in product_c_atoms.items():
                        if map_num1 < map_num2:  # Avoid checking the same pair twice
                            # Check if these atoms are bonded in product
                            bond = product_mol.GetBondBetweenAtoms(idx1, idx2)
                            if bond is not None:
                                # Find these atoms in reactants
                                atoms_found = [False, False]
                                atoms_in_same_reactant = False

                                for r_mol in reactant_mols:
                                    r_atom1_idx = None
                                    r_atom2_idx = None

                                    for atom in r_mol.GetAtoms():
                                        if atom.GetAtomMapNum() == map_num1:
                                            atoms_found[0] = True
                                            r_atom1_idx = atom.GetIdx()
                                        elif atom.GetAtomMapNum() == map_num2:
                                            atoms_found[1] = True
                                            r_atom2_idx = atom.GetIdx()

                                    # If both atoms are in this reactant, check if they're bonded
                                    if r_atom1_idx is not None and r_atom2_idx is not None:
                                        atoms_in_same_reactant = True
                                        r_bond = r_mol.GetBondBetweenAtoms(r_atom1_idx, r_atom2_idx)
                                        if r_bond is not None:
                                            # Bond exists in reactant, not a new bond
                                            break

                                # If atoms were found in reactants but not bonded, it's a new C-C bond
                                if (
                                    atoms_found[0]
                                    and atoms_found[1]
                                    and atoms_in_same_reactant == False
                                ):
                                    print(
                                        f"Detected new C-C bond formation between mapped atoms {map_num1} and {map_num2}"
                                    )
                                    new_cc_bond_found = True
                                    break

                    if new_cc_bond_found:
                        has_late_stage_cc_bond = True
                        break

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_late_stage_cc_bond
