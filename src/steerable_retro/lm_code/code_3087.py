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
    Detects if the synthesis route involves a halogen exchange (particularly between F, Cl, Br, I).
    """
    has_halogen_exchange = False

    def dfs_traverse(node):
        nonlocal has_halogen_exchange

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Check for known halogen exchange reactions
            if checker.check_reaction("Finkelstein reaction", rsmi):
                print(f"Detected Finkelstein reaction (halogen exchange): {rsmi}")
                has_halogen_exchange = True
                return

            if checker.check_reaction("Aromatic substitution of bromine by chlorine", rsmi):
                print(f"Detected aromatic halogen exchange (Br to Cl): {rsmi}")
                has_halogen_exchange = True
                return

            # For other potential halogen exchanges, analyze reactants and products
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules for proper analysis
            try:
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi.strip()]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if not all(reactant_mols) or not product_mol:
                    print(f"Warning: Could not parse all molecules in reaction: {rsmi}")
                    return

                # Check for halogen atoms in reactants and product
                halogen_patterns = {"F": "[F]", "Cl": "[Cl]", "Br": "[Br]", "I": "[I]"}

                reactant_halogens = set()
                for mol in reactant_mols:
                    for halogen, pattern in halogen_patterns.items():
                        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                            reactant_halogens.add(halogen)

                product_halogens = set()
                for halogen, pattern in halogen_patterns.items():
                    if product_mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                        product_halogens.add(halogen)

                # Check if there's a difference in halogens between reactants and product
                if reactant_halogens and product_halogens and reactant_halogens != product_halogens:
                    print(
                        f"Detected potential halogen exchange: {reactant_halogens} to {product_halogens} in {rsmi}"
                    )
                    has_halogen_exchange = True
            except Exception as e:
                print(f"Error analyzing reaction for halogen exchange: {e}")

        # Traverse children
        for child in node.get("children", []):
            if not has_halogen_exchange:  # Stop traversal if we already found a halogen exchange
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Halogen exchange strategy detected: {has_halogen_exchange}")
    return has_halogen_exchange
