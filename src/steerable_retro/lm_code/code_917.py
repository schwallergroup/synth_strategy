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
    Detects if the synthesis uses a late-stage SNAr coupling with an aniline derivative.
    Looks for a reaction where an amine attacks a carbon adjacent to electron-withdrawing groups.
    """
    found = False

    def dfs_traverse(node, depth=0):
        nonlocal found

        if found:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0-2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is an SNAr reaction using specific reaction types
                is_snar = (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                )

                if is_snar:
                    print("Found SNAr or N-arylation reaction pattern")

                    # Check for amine nucleophile in reactants
                    amine_found = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                        ):
                            print(f"Found amine nucleophile in reactant: {reactant}")
                            amine_found = True
                            break

                    # Check for aromatic halide in reactants (typical leaving group in SNAr)
                    aromatic_halide_found = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            print(f"Found aromatic halide in reactant: {reactant}")
                            aromatic_halide_found = True
                            break

                    # If we have both components of an SNAr, mark as found
                    if amine_found and aromatic_halide_found:
                        print("Confirmed late-stage SNAr coupling with amine derivative")
                        found = True
                        return

                # Alternative detection method if reaction checker doesn't identify it
                if not found:
                    amine_reactant = None
                    aromatic_halide_reactant = None

                    # Find amine and aromatic halide reactants
                    for reactant in reactants:
                        if (
                            checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                        ):
                            amine_reactant = reactant
                        if checker.check_fg("Aromatic halide", reactant):
                            aromatic_halide_reactant = reactant

                    # If we have both components, check for electron-withdrawing groups
                    if amine_reactant and aromatic_halide_reactant:
                        # Check for electron-withdrawing groups that activate SNAr
                        ewg_found = False
                        ewg_groups = [
                            "Nitro group",
                            "Nitrile",
                            "Trifluoro group",
                            "Trichloro group",
                            "Ketone",
                            "Ester",
                            "Carboxylic acid",
                            "Sulfone",
                        ]

                        for ewg in ewg_groups:
                            if checker.check_fg(ewg, aromatic_halide_reactant):
                                print(
                                    f"Found electron-withdrawing group ({ewg}) in reactant: {aromatic_halide_reactant}"
                                )
                                ewg_found = True
                                break

                        # Check if the product contains a new C-N bond (indicating coupling)
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and (ewg_found or depth <= 1):
                            # For very late-stage reactions (depth 0-1), we're more lenient
                            print(
                                "Confirmed late-stage SNAr coupling with amine derivative (alternative method)"
                            )
                            found = True
                            return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found
