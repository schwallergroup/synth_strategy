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
    This function detects if the synthesis route includes conversion of
    benzyl alcohol to benzyl chloride.
    """
    conversion_found = False

    def dfs_traverse(node, depth=0):
        nonlocal conversion_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product_part = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, we're looking for the reverse reaction:
                # Product should have benzyl chloride, reactants should have benzyl alcohol

                # Check for alcohol to chloride conversion reactions
                alcohol_to_chloride_reactions = [
                    "Alcohol to chloride_sulfonyl chloride",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                ]

                for reaction_type in alcohol_to_chloride_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"  Detected '{reaction_type}' reaction")

                        # Check if product contains benzyl chloride (in retrosynthesis)
                        product_mol = Chem.MolFromSmiles(product_part)
                        benzyl_cl_pattern = Chem.MolFromSmarts("c[CH2]Cl")
                        if product_mol and product_mol.HasSubstructMatch(benzyl_cl_pattern):
                            print(f"  Detected benzyl chloride in product: {product_part}")

                            # Check if any reactant contains benzyl alcohol
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                benzyl_oh_pattern = Chem.MolFromSmarts("c[CH2]O")
                                if reactant_mol and reactant_mol.HasSubstructMatch(
                                    benzyl_oh_pattern
                                ):
                                    print(f"  Detected benzyl alcohol in reactant: {reactant}")
                                    print(
                                        f"  Benzyl alcohol to benzyl chloride conversion detected at depth {depth}"
                                    )
                                    conversion_found = True
                                    return

                # Alternative approach: check for functional groups directly
                if not conversion_found:
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("c[CH2]Cl")
                    ):
                        print(f"  Found benzyl chloride in product: {product_part}")

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("c[CH2]O")
                            ):
                                print(f"  Found benzyl alcohol in reactant: {reactant}")
                                print(
                                    f"  Benzyl alcohol to benzyl chloride conversion detected at depth {depth}"
                                )
                                conversion_found = True
                                return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return conversion_found
