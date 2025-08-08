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
    Detects if the synthesis follows a linear strategy without convergent steps.

    A linear synthesis has sequential transformations where each reaction has only
    one significant reactant contributing to the product's core structure.
    """
    is_linear = True

    def is_simple_reagent(smiles, reaction_type=None):
        """
        Identifies if a molecule is likely a simple reagent rather than a key reactant.
        Takes into account the reaction context when available.
        """
        # Strip atom mapping before length check
        unmapped_smiles = re.sub(r":[0-9]+", "", smiles)

        # Very short SMILES are often simple reagents
        if len(unmapped_smiles) < 8:
            return True

        # Check for common reagent patterns
        common_reagents = [
            "O",
            "O=O",
            "Cl",
            "Br",
            "I",
            "F",
            "H",
            "[H]",
            "[OH-]",
            "[H+]",
            "CC(=O)O",
            "CCO",
            "CO",
            "CN",
            "CS",
            "C(=O)O",
            "C#N",
            "N",
            "S",
            "P",
            "B",
            "Si",
            "Na",
            "K",
            "Li",
            "Mg",
            "CC(=O)[O-]",
            "CC(C)O",
            "CCN",
            "CC#N",
            "C[Si](C)(C)C",
            "C[Si](C)(C)Cl",
        ]

        if unmapped_smiles in common_reagents:
            return True

        # Create molecule object
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False  # If we can't parse the SMILES, assume it's not a simple reagent

        # Count atoms excluding H
        heavy_atom_count = mol.GetNumHeavyAtoms()

        # Small molecules are often reagents
        if heavy_atom_count <= 5:
            return True

        # Check for common reagent functional groups
        if heavy_atom_count <= 10:
            # Common acid/base reagents
            if checker.check_fg("Carboxylic acid", smiles):
                return True
            if checker.check_fg("Primary alcohol", smiles) and heavy_atom_count <= 8:
                return True
            if checker.check_fg("Primary amine", smiles) and heavy_atom_count <= 8:
                return True

            # Common organometallic reagents
            if "MgBr" in unmapped_smiles or "MgCl" in unmapped_smiles or "Li" in unmapped_smiles:
                return True

            # Common protecting groups and derivatizing agents
            if "Si(C)(C)C" in unmapped_smiles:  # TMS groups
                return True
            if "C(=O)OC(C)(C)C" in unmapped_smiles:  # Boc groups
                return True
            if "S(=O)(=O)" in unmapped_smiles and heavy_atom_count <= 10:  # Sulfonyl groups
                return True

        # Check reaction context if available
        if reaction_type:
            # For certain reaction types, specific reactants are typically reagents
            if reaction_type == "Methylation" and "C" in unmapped_smiles and heavy_atom_count <= 6:
                return True
            if "Oxidation" in reaction_type and (
                "O" in unmapped_smiles or "O=O" in unmapped_smiles
            ):
                return True
            if "Reduction" in reaction_type and heavy_atom_count <= 8:
                return True
            if "Protection" in reaction_type and heavy_atom_count <= 15:
                return True

            # Special case for dimethyl acetal/ketal formation
            if (
                checker.check_reaction("Aldehyde or ketone acetalization", reaction_type)
                and "CO" in unmapped_smiles
            ):
                return True

            # Special case for Grignard reactions
            if checker.check_reaction("Grignard", reaction_type) and (
                "MgBr" in unmapped_smiles or "MgCl" in unmapped_smiles
            ):
                return True

        # Check for molecules without rings - often reagents if small enough
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 0 and heavy_atom_count <= 8:
            return True

        # If none of the above conditions are met, it's likely a significant reactant
        return False

    def get_reaction_type(node):
        """Extract reaction type from metadata if available"""
        if "metadata" in node and "ID" in node["metadata"]:
            return node["metadata"]["ID"]
        return None

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found it's not linear
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Get reaction type if available
            reaction_type = get_reaction_type(node)

            # Filter out simple reagents
            significant_reactants = [
                r for r in reactants if not is_simple_reagent(r, reaction_type)
            ]

            # If more than one significant reactant, it's convergent
            if len(significant_reactants) > 1:
                print(
                    f"Found convergent step at depth {depth} with {len(significant_reactants)} significant reactants"
                )
                print(f"Reactants: {significant_reactants}")
                print(f"Reaction type: {reaction_type}")

                # Special case for dimethyl acetal/ketal reactions
                if reaction_type and "acetal" in reaction_type.lower():
                    # These often appear convergent but are considered linear
                    print(
                        "Acetal/ketal formation detected - treating as linear despite multiple reactants"
                    )
                else:
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
