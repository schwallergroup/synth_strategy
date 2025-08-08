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
    Detects if the route involves C-N bond disconnection in the retrosynthetic direction.
    """
    cn_disconnection_found = False

    def dfs_traverse(node):
        nonlocal cn_disconnection_found

        if cn_disconnection_found:
            return  # Early return if already found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for common C-N bond forming reactions
                cn_reaction_types = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                ]

                for rxn_type in cn_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"C-N bond disconnection detected: {rxn_type}")
                        cn_disconnection_found = True
                        return

                # If no specific reaction type matched, check for C-N bond breaking using atom mapping
                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Get all C-N bonds in the product
                cn_bonds = []
                for bond in product_mol.GetBonds():
                    a1 = product_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                    a2 = product_mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                    if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or (
                        a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6
                    ):
                        # Get atom maps if available
                        a1_map = (
                            a1.GetProp("molAtomMapNumber")
                            if a1.HasProp("molAtomMapNumber")
                            else None
                        )
                        a2_map = (
                            a2.GetProp("molAtomMapNumber")
                            if a2.HasProp("molAtomMapNumber")
                            else None
                        )
                        if a1_map and a2_map:
                            cn_bonds.append((a1_map, a2_map))

                # If we have mapped C-N bonds in product, check if they exist in reactants
                if cn_bonds:
                    # Check if any C-N bond in product doesn't exist in reactants
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        for bond in cn_bonds:
                            c_map, n_map = bond
                            # Find atoms with these mapping numbers
                            c_atom = None
                            n_atom = None
                            for atom in reactant_mol.GetAtoms():
                                if atom.HasProp("molAtomMapNumber"):
                                    map_num = atom.GetProp("molAtomMapNumber")
                                    if map_num == c_map and atom.GetAtomicNum() == 6:
                                        c_atom = atom
                                    elif map_num == n_map and atom.GetAtomicNum() == 7:
                                        n_atom = atom

                            # If both mapped atoms exist in this reactant, check if they're bonded
                            if c_atom and n_atom:
                                bond_exists = reactant_mol.GetBondBetweenAtoms(
                                    c_atom.GetIdx(), n_atom.GetIdx()
                                )
                                if not bond_exists:
                                    print(f"C-N bond disconnection detected via atom mapping")
                                    cn_disconnection_found = True
                                    return

                # Check for common nitrogen-containing functional groups that might indicate C-N disconnection
                fg_checks = [
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                    "Primary amine",
                    "Secondary amine",
                    "Tertiary amine",
                    "Aniline",
                    "Urea",
                    "Carbamate",
                ]

                product_has_n_fg = False
                for fg in fg_checks:
                    if checker.check_fg(fg, product):
                        product_has_n_fg = True
                        break

                if product_has_n_fg:
                    # Check if we have multiple fragments and at least one has nitrogen
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]
                    if len(reactant_mols) >= 2:
                        n_fragments = 0
                        for mol in reactant_mols:
                            for atom in mol.GetAtoms():
                                if atom.GetAtomicNum() == 7:
                                    n_fragments += 1
                                    break

                        # If we have multiple fragments and at least one has nitrogen, likely C-N disconnection
                        if n_fragments >= 1 and len(reactant_mols) > n_fragments:
                            print("C-N bond disconnection detected based on fragment analysis")
                            cn_disconnection_found = True
                            return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_disconnection_found
