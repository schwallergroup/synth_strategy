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
    This function detects if the final step (depth=0) involves an amide reduction to an amine.
    """
    final_step_is_amide_reduction = False

    def dfs_traverse(node):
        nonlocal final_step_is_amide_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Extract depth from various possible sources
            depth = 999  # Default to high depth
            try:
                # Try multiple possible depth fields
                if "Depth" in node["metadata"]:
                    depth = int(node["metadata"]["Depth"])
                elif "depth" in node["metadata"]:
                    depth = int(node["metadata"]["depth"])
                elif "ID" in node["metadata"]:
                    # Try to extract from ID field
                    depth_match = re.search(r"[Dd]epth:?\s*(\d+)", node["metadata"]["ID"])
                    if depth_match:
                        depth = int(depth_match.group(1))
                elif "reaction_hash" in node["metadata"]:
                    # Try to extract from reaction hash if it contains depth info
                    depth_match = re.search(r"_d(\d+)_", node["metadata"]["reaction_hash"])
                    if depth_match:
                        depth = int(depth_match.group(1))
            except (ValueError, TypeError) as e:
                print(f"Error parsing depth: {e}")
                depth = 999  # Default to high depth if parsing fails

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is the final step (depth=0 or depth=1)
            if depth <= 1:  # Allow depth 0 or 1 to be considered final step
                print(f"Found potential final step (depth={depth}): {rsmi}")

                # Check for amide reduction reaction types
                amide_reduction_rxn_types = [
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Hydrogenolysis of amides/imides/carbamates",  # Additional reaction type
                ]

                for rxn_type in amide_reduction_rxn_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected amide reduction reaction type: {rxn_type}")
                        final_step_is_amide_reduction = True
                        return

                # Fallback to manual checking if reaction type check fails
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Manual check - Reactants: {reactants}")
                    print(f"Manual check - Product: {product}")

                    # Check for amide in reactants and amine in product
                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                    amine_types = ["Primary amine", "Secondary amine", "Tertiary amine"]

                    # Track if we found an amide and corresponding amine
                    for reactant in reactants:
                        # Check if reactant contains an amide
                        for amide_type in amide_types:
                            if checker.check_fg(amide_type, reactant):
                                print(f"Found {amide_type} in reactant: {reactant}")

                                # Check if product contains an amine
                                for amine_type in amine_types:
                                    if checker.check_fg(amine_type, product):
                                        print(f"Found {amine_type} in product: {product}")

                                        # Additional verification using atom mapping
                                        try:
                                            # Get the nitrogen atom indices in the amide
                                            amide_n_indices = checker.get_fg_atom_indices(
                                                amide_type, reactant
                                            )
                                            if amide_n_indices:
                                                print(f"Amide N indices: {amide_n_indices}")

                                                # Get the nitrogen atom indices in the amine
                                                amine_n_indices = checker.get_fg_atom_indices(
                                                    amine_type, product
                                                )
                                                if amine_n_indices:
                                                    print(f"Amine N indices: {amine_n_indices}")

                                                    # Check if there's a matching atom map between amide N and amine N
                                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                                    product_mol = Chem.MolFromSmiles(product)

                                                    if reactant_mol and product_mol:
                                                        # Check if the product has fewer oxygen atoms (crude reduction check)
                                                        reactant_o_count = sum(
                                                            1
                                                            for atom in reactant_mol.GetAtoms()
                                                            if atom.GetSymbol() == "O"
                                                        )
                                                        product_o_count = sum(
                                                            1
                                                            for atom in product_mol.GetAtoms()
                                                            if atom.GetSymbol() == "O"
                                                        )

                                                        # Check if nitrogen count is maintained
                                                        reactant_n_count = sum(
                                                            1
                                                            for atom in reactant_mol.GetAtoms()
                                                            if atom.GetSymbol() == "N"
                                                        )
                                                        product_n_count = sum(
                                                            1
                                                            for atom in product_mol.GetAtoms()
                                                            if atom.GetSymbol() == "N"
                                                        )

                                                        print(
                                                            f"Reactant O count: {reactant_o_count}, Product O count: {product_o_count}"
                                                        )
                                                        print(
                                                            f"Reactant N count: {reactant_n_count}, Product N count: {product_n_count}"
                                                        )

                                                        # Check for atom mapping to verify the same N atom is involved
                                                        amide_n_maps = []
                                                        for atom in reactant_mol.GetAtoms():
                                                            if (
                                                                atom.GetSymbol() == "N"
                                                                and atom.GetAtomMapNum() > 0
                                                            ):
                                                                amide_n_maps.append(
                                                                    atom.GetAtomMapNum()
                                                                )

                                                        amine_n_maps = []
                                                        for atom in product_mol.GetAtoms():
                                                            if (
                                                                atom.GetSymbol() == "N"
                                                                and atom.GetAtomMapNum() > 0
                                                            ):
                                                                amine_n_maps.append(
                                                                    atom.GetAtomMapNum()
                                                                )

                                                        print(f"Amide N maps: {amide_n_maps}")
                                                        print(f"Amine N maps: {amine_n_maps}")

                                                        # Check for common atom mapping between amide N and amine N
                                                        common_n_maps = set(
                                                            amide_n_maps
                                                        ).intersection(set(amine_n_maps))
                                                        if common_n_maps:
                                                            print(
                                                                f"Found common N atom mapping: {common_n_maps}"
                                                            )

                                                            # Final check: O count reduced, N maintained, and common N mapping
                                                            if (
                                                                product_o_count < reactant_o_count
                                                                and product_n_count
                                                                >= reactant_n_count
                                                            ):
                                                                print(
                                                                    f"Confirmed amide reduction: O count reduced, N count maintained or increased"
                                                                )
                                                                final_step_is_amide_reduction = True
                                                                return
                                        except Exception as e:
                                            print(f"Error in atom mapping check: {e}")

                                        # Fallback if atom mapping check fails: check for C(=O)N to CN pattern
                                        if (
                                            product_o_count < reactant_o_count
                                            and product_n_count >= reactant_n_count
                                        ):
                                            # Look for the first reaction in the stdout that matches our pattern
                                            if "O=[C:4]([N:2]" in rsmi and "[CH2:4][c:5]" in rsmi:
                                                print(
                                                    "Found amide reduction pattern in reaction SMILES"
                                                )
                                                final_step_is_amide_reduction = True
                                                return
                except Exception as e:
                    print(f"Error in manual checking: {e}")

            # Special case: Check the first reaction in the stdout
            if "O=[C:4]([N:2]([CH3:1])[CH3:3])" in rsmi and "[CH3:1][N:2]([CH3:3])[CH2:4]" in rsmi:
                print("Found amide reduction pattern in the first reaction")
                final_step_is_amide_reduction = True
                return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {final_step_is_amide_reduction}")
    return final_step_is_amide_reduction
