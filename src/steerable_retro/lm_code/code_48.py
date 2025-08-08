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
    This function detects late-stage SNAr coupling of complex fragments.
    Specifically looking for C-O or C-N bond formation via SNAr in the final steps.
    """
    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a late-stage reaction (depth 0, 1, or 2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for SNAr reactions using the checker functions
                if (
                    checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                ):

                    print(f"SNAr coupling reaction detected at depth {depth}")
                    snar_detected = True
                    return

                # If no specific reaction type matched, check for the pattern manually
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: aromatic halide + nucleophile -> new C-O or C-N bond
                aromatic_halide_found = False
                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        print(f"Found aromatic halide in reactant: {reactant}")
                        aromatic_halide_found = True

                        # Check if other reactants contain nucleophiles (O or N containing)
                        for other_reactant in reactants:
                            if (
                                checker.check_fg("Phenol", other_reactant)
                                or checker.check_fg("Primary alcohol", other_reactant)
                                or checker.check_fg("Secondary alcohol", other_reactant)
                                or checker.check_fg("Primary amine", other_reactant)
                                or checker.check_fg("Secondary amine", other_reactant)
                                or checker.check_fg("Tertiary amine", other_reactant)
                                or checker.check_fg("Aniline", other_reactant)
                            ):
                                print(f"Found nucleophile in reactant: {other_reactant}")

                                # Check if product has new C-O or C-N bond
                                try:
                                    # Look for ether or amine formation in product
                                    if (
                                        checker.check_fg("Ether", product)
                                        or checker.check_fg("Aniline", product)
                                        or checker.check_fg("Primary amine", product)
                                        or checker.check_fg("Secondary amine", product)
                                        or checker.check_fg("Tertiary amine", product)
                                    ):
                                        print(f"Found potential SNAr product: {product}")
                                        snar_detected = True
                                        return
                                except:
                                    pass

                # Special case: Check if a single reactant contains both aromatic halide and nucleophile
                # This handles intramolecular SNAr
                if not aromatic_halide_found:
                    for reactant in reactants:
                        if (
                            "Cl" in reactant
                            or "Br" in reactant
                            or "I" in reactant
                            or "F" in reactant
                        ):
                            # Check if this is an aromatic system with a halide
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                for atom in reactant_mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() in ["Cl", "Br", "I", "F"]
                                        and atom.GetIsAromatic()
                                    ):
                                        print(
                                            f"Found aromatic halide in complex reactant: {reactant}"
                                        )

                                        # Check if product has new C-O or C-N bond where halide was
                                        if (
                                            checker.check_fg("Ether", product)
                                            or checker.check_fg("Aniline", product)
                                            or checker.check_fg("Primary amine", product)
                                            or checker.check_fg("Secondary amine", product)
                                            or checker.check_fg("Tertiary amine", product)
                                        ):
                                            print(f"Found potential SNAr product: {product}")
                                            snar_detected = True
                                            return

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return snar_detected
