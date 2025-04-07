#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis is convergent with exactly two fragments combined in the final step.
    """
    print("Starting two_fragment_convergent_synthesis analysis")
    is_two_fragment_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_two_fragment_convergent

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            print(f"Found reaction node at depth {depth}")

            # The final reaction step is typically at depth 1 (right below the root product)
            if depth == 1:
                print("Analyzing final reaction step")
                try:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Final reaction SMILES: {rsmi}")

                    # Extract reactants and product
                    reactants_part = rsmi.split(">")[0]
                    reactants_smiles = reactants_part.split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Found {len(reactants_smiles)} reactants")

                    # Check if there are exactly two reactants
                    if len(reactants_smiles) == 2:
                        # Create RDKit molecules for analysis
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        # Ensure both reactants are substantial fragments (at least 6 atoms)
                        if (
                            all(reactant_mols)
                            and product_mol
                            and all(len(mol.GetAtoms()) > 5 for mol in reactant_mols)
                        ):
                            print(
                                f"Found two substantial fragments with {[len(mol.GetAtoms()) for mol in reactant_mols]} atoms"
                            )

                            # Expanded list of coupling and fragment-combining reactions
                            coupling_reactions = [
                                "Suzuki coupling with boronic acids",
                                "Suzuki coupling with boronic esters",
                                "Suzuki coupling with boronic acids OTf",
                                "Suzuki coupling with boronic esters OTf",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                                "Negishi coupling",
                                "Stille reaction_aryl",
                                "Stille reaction_vinyl",
                                "Stille reaction_benzyl",
                                "Stille reaction_allyl",
                                "Sonogashira alkyne_aryl halide",
                                "Sonogashira acetylene_aryl halide",
                                "Heck terminal vinyl",
                                "Heck_terminal_vinyl",
                                "Heck_non-terminal_vinyl",
                                "Wittig reaction with triphenylphosphorane",
                                "Wittig",
                                "Diels-Alder",
                                "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                                "Huisgen_Cu-catalyzed_1,4-subst",
                                "Huisgen_Ru-catalyzed_1,5_subst",
                                "Ullmann condensation",
                                "Ullmann-Goldberg Substitution amine",
                                "Ullmann-Goldberg Substitution thiol",
                                "Ullmann-Goldberg Substitution aryl alcohol",
                                "Hiyama-Denmark Coupling",
                                "Kumada cross-coupling",
                                "Aryllithium cross-coupling",
                                "decarboxylative_coupling",
                            ]

                            # Check if this is a known coupling reaction
                            for rxn in coupling_reactions:
                                if checker.check_reaction(rxn, rsmi):
                                    print(f"Found coupling reaction: {rxn}")
                                    is_two_fragment_convergent = True
                                    return

                            # If not a known coupling reaction, check if fragments are combined
                            # Ensure reactants are separate fragments and product is a single fragment
                            if (
                                all(len(Chem.GetMolFrags(mol)) == 1 for mol in reactant_mols)
                                and len(Chem.GetMolFrags(product_mol)) == 1
                            ):
                                print(
                                    "Reactants are separate fragments and product is a single fragment"
                                )

                                # Additional check for amide formation, ether formation, or other common linkages
                                common_linkage_reactions = [
                                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                    "Esterification of Carboxylic Acids",
                                    "Williamson Ether Synthesis",
                                    "Williamson_ether",
                                    "Mitsunobu esterification",
                                    "Mitsunobu aryl ether",
                                    "Urea synthesis via isocyanate and primary amine",
                                    "Urea synthesis via isocyanate and secondary amine",
                                    "Schotten-Baumann_amide",
                                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                    "Carboxylic acid with primary amine to amide",
                                    "Ester with primary amine to amide",
                                ]

                                for rxn in common_linkage_reactions:
                                    if checker.check_reaction(rxn, rsmi):
                                        print(f"Found fragment-linking reaction: {rxn}")
                                        is_two_fragment_convergent = True
                                        return

                                # Generic check for new bond formation between fragments
                                # If we got here and the fragments are combined, it's likely a convergent synthesis
                                print(
                                    "Found evidence of fragment combination - marking as convergent synthesis"
                                )
                                is_two_fragment_convergent = True
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final result: is_two_fragment_convergent = {is_two_fragment_convergent}")
    return is_two_fragment_convergent
