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


def main(route):
    """
    Detects if the synthetic route involves bromination followed by nucleophilic substitution,
    which is a common strategy for introducing new substituents.
    """
    # Track bromination and substitution reactions
    bromination_reactions = []
    substitution_reactions = []
    reaction_depths = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Store reaction depth
            reaction_id = node["metadata"].get("ID", str(depth))
            reaction_depths[reaction_id] = depth

            try:
                # Check for bromination: addition of Br
                product_mol = Chem.MolFromSmiles(product)

                has_bromo_in_product = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]-[Br]")
                )

                has_bromo_in_reactants = False
                for reactant in reactants:
                    if "Br" in reactant:  # Simple check for bromine-containing reagent
                        has_bromo_in_reactants = True
                        break

                if has_bromo_in_product and has_bromo_in_reactants:
                    bromination_reactions.append(reaction_id)

                # Check for nucleophilic substitution: Br replaced by N, O, S
                has_bromo_in_reactants_mol = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[Br]")
                    ):
                        has_bromo_in_reactants_mol = True
                        break

                # Check if product has new C-N, C-O or C-S bond
                if has_bromo_in_reactants_mol and product_mol:
                    if (
                        product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#7]"))
                        or product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#8]"))
                        or product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#16]"))
                    ):
                        substitution_reactions.append(reaction_id)
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's a bromination followed by substitution
    for bromination_id in bromination_reactions:
        for substitution_id in substitution_reactions:
            # If substitution depth is less than bromination depth, it comes after in synthesis
            if reaction_depths[substitution_id] < reaction_depths[bromination_id]:
                print(
                    f"Found bromination-substitution sequence: bromination at depth {reaction_depths[bromination_id]}, substitution at depth {reaction_depths[substitution_id]}"
                )
                return True

    return False
