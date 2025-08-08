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
    Detects a strategy involving N-demethylation as one of the final steps
    """
    found_demethylation = False
    demethylation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_demethylation, demethylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                prod_mol = Chem.MolFromSmiles(product)

                # Check for N-demethylation pattern
                for reactant in reactants:
                    react_mol = Chem.MolFromSmiles(reactant)

                    if react_mol and prod_mol:
                        # Since we're traversing retrosynthetically, in a demethylation:
                        # - The product (in rsmi) should have a tertiary amine
                        # - The reactant (in rsmi) should have a secondary amine
                        has_secondary_amine = checker.check_fg("Secondary amine", reactant)
                        has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                        if has_secondary_amine and has_tertiary_amine:
                            # Check if this is a methylation reaction (which means demethylation in retrosynthesis)
                            is_demethylation = False

                            # Check for specific methylation reaction types
                            methylation_reactions = [
                                "N-methylation",
                                "Eschweiler-Clarke Secondary Amine Methylation",
                                "Eschweiler-Clarke Primary Amine Methylation",
                                "Reductive methylation of primary amine with formaldehyde",
                                "Methylation with MeI_primary",
                                "Methylation with MeI_secondary",
                                "Methylation with MeI_tertiary",
                                "DMS Amine methylation",
                                "Hydrogenolysis of tertiary amines",
                                "Parnes methylation",
                            ]

                            for reaction_type in methylation_reactions:
                                if checker.check_reaction(reaction_type, rsmi):
                                    is_demethylation = True
                                    print(f"Detected {reaction_type} reaction")
                                    break

                            # Additional structural check if reaction type check fails
                            if not is_demethylation:
                                # Check for methyl group attached to nitrogen in product but not in reactant
                                methyl_n_patt = Chem.MolFromSmarts("[CH3]-[N]")
                                if prod_mol.HasSubstructMatch(
                                    methyl_n_patt
                                ) and not react_mol.HasSubstructMatch(methyl_n_patt):
                                    # Verify this is likely a demethylation by checking atom counts
                                    if prod_mol.GetNumAtoms() > react_mol.GetNumAtoms():
                                        is_demethylation = True
                                        print("Detected structural pattern for N-demethylation")

                            if is_demethylation:
                                print(f"Found N-demethylation at depth {depth}, rsmi: {rsmi}")
                                found_demethylation = True
                                demethylation_depth = min(demethylation_depth, depth)
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if demethylation occurs at a low depth (late in synthesis)
    # Consider depths 0, 1, or 2 as late-stage
    if found_demethylation and demethylation_depth <= 2:
        print(f"Found late-stage demethylation strategy at depth {demethylation_depth}")
        return True

    return False
