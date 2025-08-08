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
    Detects if the synthesis route uses reductive amination in the late stage (low depth).
    Looks for C=O + HNR2 → C-NR2 transformation.
    """
    found = False

    def dfs_traverse(node, depth=0):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a low depth reaction (late stage)
            current_depth = depth
            if "ID" in node.get("metadata", {}):
                depth_str = node["metadata"]["ID"]
                if "Depth:" in depth_str:
                    try:
                        current_depth = int(depth_str.split("Depth:")[1].strip().split()[0])
                    except:
                        print(
                            f"Could not parse depth from {depth_str}, using traversal depth {depth}"
                        )

            # Consider depth 0, 1, and 2 as late stage
            if current_depth <= 2:
                print(f"Checking reaction at depth {current_depth}: {rsmi}")

                # Use checker to detect reductive amination reactions
                is_reductive_amination = (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                )

                if is_reductive_amination:
                    print(f"Found late-stage reductive amination at depth {current_depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    found = True
                    return

                # Fallback method if checker doesn't detect it
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carbonyl groups in reactants
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                has_carbonyl = has_aldehyde or has_ketone

                # Check for amines in reactants
                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)
                has_amine = has_primary_amine or has_secondary_amine

                if has_carbonyl and has_amine:
                    # Check if product has a new C-N bond that wasn't in reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if product_mol is not None:
                        # Check for C-N bond pattern in product
                        cn_bond_pattern = Chem.MolFromSmarts("[C]-[N]")

                        # Count C-N bonds in reactants and product
                        reactant_cn_count = sum(
                            len(mol.GetSubstructMatches(cn_bond_pattern))
                            for mol in reactant_mols
                            if mol is not None
                        )

                        product_cn_count = len(product_mol.GetSubstructMatches(cn_bond_pattern))

                        # Check if there are more C-N bonds in the product than in the reactants
                        if product_cn_count > reactant_cn_count:
                            print(
                                f"Found late-stage reductive amination at depth {current_depth} (fallback method)"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            print(f"Carbonyl: {has_carbonyl}, Amine: {has_amine}")
                            print(
                                f"C-N bonds in reactants: {reactant_cn_count}, in product: {product_cn_count}"
                            )
                            found = True
                            return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late-stage reductive amination found: {found}")
    return found
