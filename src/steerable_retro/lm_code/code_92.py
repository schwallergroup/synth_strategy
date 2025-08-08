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
    This function detects late-stage halogen manipulation, particularly
    dehalogenation or halogenation in the final steps of synthesis.
    """
    has_late_halogen_manipulation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_halogen_manipulation

        if node["type"] == "reaction" and depth <= 2:  # Focus on late-stage reactions (low depth)
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a dehalogenation reaction using reaction checkers
                dehalogenation_reactions = ["Aromatic dehalogenation", "Dehalogenation"]

                # Check for halogenation reactions
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Fluorination",
                    "Chlorination",
                    "Bromination",
                    "Iodination",
                    "Aromatic substitution of bromine by chlorine",
                    "Halodeboronation of boronic acids",
                    "Halodeboronation of boronic esters",
                ]

                # Check for known reaction types first
                for rxn_type in dehalogenation_reactions + halogenation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_late_halogen_manipulation = True
                        print(f"Detected late-stage {rxn_type} reaction at depth {depth}")
                        return

                # If not a known reaction type, check for halogen manipulation manually
                if product_smiles and Chem.MolFromSmiles(product_smiles):
                    # Define all halogen functional groups to check
                    halogen_fgs = [
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Alkenyl halide",
                        "Haloalkyne",
                    ]

                    # Count halogens in product
                    product_halogen_count = 0
                    for fg in halogen_fgs:
                        if checker.check_fg(fg, product_smiles):
                            product_halogen_count += len(
                                checker.get_fg_atom_indices(fg, product_smiles)
                            )

                    # Count halogens in reactants
                    reactants_halogen_count = 0
                    for reactant in reactants_smiles:
                        if not reactant or not Chem.MolFromSmiles(reactant):
                            continue
                        for fg in halogen_fgs:
                            if checker.check_fg(fg, reactant):
                                reactants_halogen_count += len(
                                    checker.get_fg_atom_indices(fg, reactant)
                                )

                    # Check if number of halogens changed
                    if product_halogen_count != reactants_halogen_count:
                        has_late_halogen_manipulation = True
                        if product_halogen_count > reactants_halogen_count:
                            print(f"Detected late-stage halogen addition at depth {depth}")
                        else:
                            print(f"Detected late-stage halogen removal at depth {depth}")
                        return

                    # Check for specific halogen types (F, Cl, Br, I)
                    specific_halogens = [
                        "Aromatic fluoride",
                        "Aromatic chloride",
                        "Aromatic bromide",
                        "Aromatic iodide",
                    ]

                    # Check if halogen type changed (exchange)
                    for fg in specific_halogens:
                        reactant_count = 0
                        for reactant in reactants_smiles:
                            if not reactant or not Chem.MolFromSmiles(reactant):
                                continue
                            if checker.check_fg(fg, reactant):
                                reactant_count += len(checker.get_fg_atom_indices(fg, reactant))

                        product_count = 0
                        if checker.check_fg(fg, product_smiles):
                            product_count = len(checker.get_fg_atom_indices(fg, product_smiles))

                        if reactant_count != product_count:
                            has_late_halogen_manipulation = True
                            print(f"Detected late-stage {fg} manipulation at depth {depth}")
                            return

                    # Additional check for chlorination specifically
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    for atom in product_mol.GetAtoms():
                        if (
                            atom.GetSymbol() == "Cl"
                            and atom.GetIdx() == int(atom.GetProp("molAtomMapNumber")) - 1
                        ):
                            # This is an atom-mapped chlorine, check if it was added
                            for reactant in reactants_smiles:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if not reactant_mol:
                                    continue

                                # Look for the same atom map number in reactants
                                atom_map = atom.GetProp("molAtomMapNumber")
                                found_in_reactant = False
                                for r_atom in reactant_mol.GetAtoms():
                                    if (
                                        r_atom.HasProp("molAtomMapNumber")
                                        and r_atom.GetProp("molAtomMapNumber") == atom_map
                                    ):
                                        if r_atom.GetSymbol() != "Cl":
                                            # Atom was not Cl in reactant but is Cl in product
                                            has_late_halogen_manipulation = True
                                            print(
                                                f"Detected late-stage chlorination at depth {depth}"
                                            )
                                            return
                                        found_in_reactant = True
                                        break

                                if not found_in_reactant:
                                    # New chlorine atom was added
                                    has_late_halogen_manipulation = True
                                    print(f"Detected late-stage chlorine addition at depth {depth}")
                                    return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_halogen_manipulation
