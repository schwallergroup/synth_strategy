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
    Detects if the synthesis route involves formation of a tetrasubstituted alkene
    (a C=C bond where both carbons are connected to two non-hydrogen atoms).
    """
    has_tetrasubstituted_alkene_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tetrasubstituted_alkene_formation

        if node["type"] == "reaction" and not has_tetrasubstituted_alkene_formation:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth} - Analyzing reaction: {rsmi}")

                # Check for reactions known to form tetrasubstituted alkenes
                tetra_alkene_forming_reactions = [
                    "Wittig",
                    "Wittig with Phosphonium",
                    "Julia Olefination",
                    "Aldol condensation",
                    "Heck",
                    "Olefination of ketones with Grignard reagents",
                    "Olefination of aldehydes with Grignard reagents",
                ]

                for rxn_type in tetra_alkene_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(
                            f"Detected potential tetrasubstituted alkene-forming reaction: {rxn_type}"
                        )

                        # Now verify that a tetrasubstituted alkene was actually formed
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]

                        if product_mol and all(m for m in reactant_mols):
                            # More general pattern for tetrasubstituted alkene
                            tetra_alkene_pattern = Chem.MolFromSmarts(
                                "[#6](-[!#1])(-[!#1])=[#6](-[!#1])(-[!#1])"
                            )

                            # Check if product has tetrasubstituted alkene
                            product_matches = product_mol.GetSubstructMatches(tetra_alkene_pattern)
                            if product_matches:
                                print(f"Product contains tetrasubstituted alkene: {product_smiles}")

                                # Check if any reactant already has the same tetrasubstituted alkene
                                reactant_has_tetra_alkene = False
                                for r_mol in reactant_mols:
                                    if r_mol.GetSubstructMatches(tetra_alkene_pattern):
                                        reactant_has_tetra_alkene = True
                                        print(f"Reactant already contains tetrasubstituted alkene")
                                        break

                                if not reactant_has_tetra_alkene:
                                    has_tetrasubstituted_alkene_formation = True
                                    print(
                                        f"Confirmed tetrasubstituted alkene formation in reaction: {rsmi}"
                                    )
                                    return

                # If no specific reaction type was detected, try a more general approach
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]

                if product_mol and all(m for m in reactant_mols):
                    # Pattern for tetrasubstituted alkene
                    tetra_alkene_pattern = Chem.MolFromSmarts(
                        "[#6](-[!#1])(-[!#1])=[#6](-[!#1])(-[!#1])"
                    )

                    # Check if product has tetrasubstituted alkene
                    product_matches = product_mol.GetSubstructMatches(tetra_alkene_pattern)
                    if product_matches:
                        print(f"Product contains tetrasubstituted alkene: {product_smiles}")

                        # Check if reactants don't have the same tetrasubstituted alkene
                        reactant_has_tetra_alkene = False
                        for r_mol in reactant_mols:
                            if r_mol.GetSubstructMatches(tetra_alkene_pattern):
                                reactant_has_tetra_alkene = True
                                print(f"Reactant already contains tetrasubstituted alkene")
                                break

                        if not reactant_has_tetra_alkene:
                            # Check for functional groups that suggest tetrasubstituted alkene formation
                            has_relevant_fg = False

                            # Check for carbonyl groups in reactants
                            for r_smi in reactants_smiles:
                                if checker.check_fg("Aldehyde", r_smi) or checker.check_fg(
                                    "Ketone", r_smi
                                ):
                                    has_relevant_fg = True
                                    print(f"Reactant contains carbonyl group: {r_smi}")
                                    break

                            # Check for phosphorus ylides or other olefination reagents
                            for r_smi in reactants_smiles:
                                if "P" in r_smi:
                                    has_relevant_fg = True
                                    print(f"Reactant contains phosphorus: {r_smi}")
                                    break

                            if has_relevant_fg:
                                has_tetrasubstituted_alkene_formation = True
                                print(
                                    f"Confirmed tetrasubstituted alkene formation based on functional groups"
                                )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if not has_tetrasubstituted_alkene_formation:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_tetrasubstituted_alkene_formation
