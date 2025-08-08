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
    This function detects if the synthesis follows a convergent strategy with
    multiple fragments being combined in late-stage reactions.
    """
    # Track reactions and their properties
    reaction_data = []

    # List of coupling reactions commonly used in convergent synthesis
    coupling_reactions = [
        "Suzuki",
        "Negishi",
        "Heck",
        "Sonogashira",
        "Stille",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg",
        "Kumada",
    ]

    # List of significant ring structures
    significant_rings = [
        "benzene",
        "pyridine",
        "naphthalene",
        "indole",
        "quinoline",
        "isoquinoline",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
    ]

    def is_significant_fragment(smiles):
        """Determine if a molecule is a significant fragment"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check size - consider molecules with more than 7 heavy atoms significant
        if mol.GetNumHeavyAtoms() > 7:
            return True

        # Check for significant ring structures
        for ring in significant_rings:
            if checker.check_ring(ring, smiles):
                return True

        return False

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count significant fragments
                significant_reactants = [r for r in reactants if is_significant_fragment(r)]

                # Check if this is a coupling reaction
                is_coupling = any(checker.check_reaction(rxn, rsmi) for rxn in coupling_reactions)

                # Store reaction data for analysis
                reaction_data.append(
                    {
                        "depth": depth,
                        "reactants": reactants,
                        "significant_reactants": significant_reactants,
                        "is_coupling": is_coupling,
                        "rsmi": rsmi,
                    }
                )

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze reaction data to determine if synthesis is convergent
    convergent = False

    # Look for late-stage reactions (depth <= 3) that combine multiple significant fragments
    for reaction in reaction_data:
        if reaction["depth"] <= 3:
            # Check if this reaction combines multiple significant fragments
            if len(reaction["significant_reactants"]) >= 2:
                print(
                    f"Found convergent reaction at depth {reaction['depth']} with {len(reaction['significant_reactants'])} significant fragments"
                )

                # If it's a known coupling reaction, that's strong evidence of convergent synthesis
                if reaction["is_coupling"]:
                    print(f"  This is a coupling reaction: {reaction['rsmi']}")
                    convergent = True
                    break

                # Even if not a known coupling, combining multiple significant fragments is evidence of convergent synthesis
                print(f"  Multiple significant fragments combined: {reaction['rsmi']}")
                convergent = True
                break

    print(f"Convergent synthesis detected: {convergent}")
    return convergent
