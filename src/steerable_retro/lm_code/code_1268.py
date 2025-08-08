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
    Detects if the synthetic route employs sequential methylation reactions
    to build molecular complexity.
    """
    methylation_reactions = []

    def dfs_traverse(node, depth):
        nonlocal methylation_reactions

        # Add depth to node for reference
        node["depth"] = depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for methylation reactions using the checker function
                methylation_types = [
                    "Methylation",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "Methylation with MeI_aryl",
                    "Methylation with MeI_SH",
                    "Methylation with DMS",
                    "Methylation with DMC",
                    "DMS COOH methylation",
                    "DMS Amine methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "Parnes methylation",
                    "N-methylation",
                    "S-methylation",
                    "O-methylation",
                    "C-methylation",
                ]

                for rxn_type in methylation_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        methylation_reactions.append((depth, rsmi, rxn_type))
                        print(f"Found methylation reaction ({rxn_type}) at depth {depth}: {rsmi}")
                        break

                # If no specific methylation reaction type was found, check structurally
                if not any(
                    checker.check_reaction(rxn_type, rsmi) for rxn_type in methylation_types
                ):
                    try:
                        # Extract reactants and product
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant contains a methyl group (CH3)
                        methyl_reactant = any("C" in reactant for reactant in reactants)

                        # Check if the reaction involves adding a methyl group
                        if methyl_reactant and product:
                            # Simple heuristic: if a reactant has "C" and the product has more carbon atoms
                            prod_mol = Chem.MolFromSmiles(product)
                            if prod_mol:
                                # Count carbon atoms in product
                                carbon_count_prod = sum(
                                    1 for atom in prod_mol.GetAtoms() if atom.GetSymbol() == "C"
                                )

                                # Check if any reactant has fewer carbon atoms
                                for reactant in reactants:
                                    react_mol = Chem.MolFromSmiles(reactant)
                                    if react_mol:
                                        carbon_count_react = sum(
                                            1
                                            for atom in react_mol.GetAtoms()
                                            if atom.GetSymbol() == "C"
                                        )
                                        if carbon_count_prod > carbon_count_react:
                                            # This might be a methylation reaction
                                            methylation_reactions.append(
                                                (depth, rsmi, "Structural methylation")
                                            )
                                            print(
                                                f"Found structural methylation at depth {depth}: {rsmi}"
                                            )
                                            break
                    except Exception as e:
                        print(f"Error analyzing reaction structure: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start DFS traversal with depth 0
    dfs_traverse(route, 0)

    # Sort by depth to check if they're sequential
    methylation_reactions.sort(key=lambda x: x[0])

    print(f"Total methylation reactions: {len(methylation_reactions)}")

    # Check if there are at least 2 methylation reactions
    if len(methylation_reactions) >= 2:
        # Check if they're at different depths (sequential in synthesis)
        for i in range(len(methylation_reactions) - 1):
            # In retrosynthesis, different depths indicate sequential steps
            if methylation_reactions[i][0] != methylation_reactions[i + 1][0]:
                print(
                    f"Found sequential methylation strategy: {methylation_reactions[i][2]} at depth {methylation_reactions[i][0]} followed by {methylation_reactions[i+1][2]} at depth {methylation_reactions[i+1][0]}"
                )
                return True

    return False
