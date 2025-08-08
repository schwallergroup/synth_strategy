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
    Detects if the synthetic route uses a convergent approach where two complex fragments
    are joined in the final step, rather than a purely linear synthesis.
    """
    convergent_strategy = False

    def is_complex_fragment(smiles):
        """Determine if a molecule is a complex fragment based on structure"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Complex fragments should have sufficient size
        if mol.GetNumHeavyAtoms() < 8:
            return False

        # Complex fragments typically contain rings
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 0:
            return False

        return True

    def is_coupling_reaction(rxn_smiles):
        """Check if the reaction is a coupling reaction"""
        # Check common coupling reactions
        coupling_reactions = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Buchwald-Hartwig",
            "N-arylation",
            "Sonogashira",
            "Heck",
            "Negishi",
            "Stille",
            "Amide formation",
            "Esterification",
            "Williamson Ether Synthesis",
        ]

        for rxn_type in coupling_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"Found coupling reaction: {rxn_type}")
                return True

        # Check for bond formation between two complex fragments
        try:
            rsmi = rxn_smiles
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # If we have multiple reactants and they're being joined
            if "." in reactants_part and not "." in product_part:
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_strategy

        # Check reactions at late stages (final or penultimate steps)
        if node["type"] == "reaction" and depth <= 1:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")

                    # Check if we have at least 2 complex reactants
                    complex_reactant_count = 0
                    for reactant in reactants:
                        if is_complex_fragment(reactant):
                            complex_reactant_count += 1

                    # Convergent synthesis requires at least 2 complex fragments
                    # and a coupling reaction
                    if complex_reactant_count >= 2 and is_coupling_reaction(rsmi):
                        print(f"Found convergent fragment coupling strategy at depth {depth}")
                        print(f"Reaction: {rsmi}")
                        print(f"Complex reactants: {complex_reactant_count}")
                        convergent_strategy = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_strategy
