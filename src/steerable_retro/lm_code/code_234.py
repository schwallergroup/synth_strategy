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
    This function detects if a stereocenter is preserved throughout the synthesis.
    It checks if stereocenters in the final product are maintained from starting materials
    through all intermediate steps.
    """
    # First check if the final product has stereocenters
    final_product = route["smiles"]
    final_mol = Chem.MolFromSmiles(final_product)
    if not final_mol:
        print(f"Could not parse final product: {final_product}")
        return False

    final_stereocenters = Chem.FindMolChiralCenters(final_mol, includeUnassigned=False)
    if not final_stereocenters:
        print(f"Final product has no stereocenters: {final_product}")
        return False

    print(f"Final product has {len(final_stereocenters)} stereocenters: {final_product}")

    # Track stereocenters through the synthesis
    preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal preserved

        if not preserved:  # Early termination if preservation is already broken
            return

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                print(f"Molecule at depth {depth}: {smiles}")
                print(f"  Stereocenters: {stereocenters}")

                # Store depth for reference
                node["depth"] = depth

                # For starting materials, we don't check preservation (they're the source)
                if node.get("in_stock", False):
                    print(f"  Starting material, no need to check preservation")
                    return
            else:
                print(f"Could not parse molecule at depth {depth}: {smiles}")

        elif node["type"] == "reaction":
            # Check if this reaction might affect stereochemistry
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction at depth {depth}: {rsmi}")

                # Check for reactions known to potentially affect stereochemistry
                stereo_affecting_reactions = [
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of aldehydes and ketones to alcohols",
                    "Oxidation of secondary alcohol to ketone",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Diels-Alder",
                    "Aldol condensation",
                ]

                for rxn_type in stereo_affecting_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"  Warning: Reaction {rxn_type} may affect stereochemistry")

                        # Check if product has expected stereocenters
                        product = rsmi.split(">")[-1]
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            product_stereocenters = Chem.FindMolChiralCenters(
                                product_mol, includeUnassigned=False
                            )

                            # Check reactants
                            reactants = rsmi.split(">")[0].split(".")
                            reactant_stereocenters = []
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    reactant_stereocenters.extend(
                                        Chem.FindMolChiralCenters(
                                            reactant_mol, includeUnassigned=False
                                        )
                                    )

                            # If reactants had stereocenters but product has fewer, flag as not preserved
                            if len(reactant_stereocenters) > 0 and len(product_stereocenters) < len(
                                reactant_stereocenters
                            ):
                                print(
                                    f"  Stereocenter count decreased: {len(reactant_stereocenters)} → {len(product_stereocenters)}"
                                )
                                preserved = False
                                return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all stereocenters were preserved
    if preserved:
        print("All stereocenters appear to be preserved throughout the synthesis")
    else:
        print("Some stereocenters were not preserved during the synthesis")

    return preserved
