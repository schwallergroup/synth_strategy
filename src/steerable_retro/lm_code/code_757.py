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
    This function detects if the synthetic route involves a late-stage deprotection step.
    Late-stage is defined as the final step or steps close to the final product (depth <= 2).
    """
    late_stage_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_deprotection_found

        # Consider depth <= 2 as late-stage
        is_late_stage = depth <= 2

        if is_late_stage and node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for various deprotection reactions
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                or checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
                or checker.check_reaction("Hydroxyl benzyl deprotection", rsmi)
                or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                or checker.check_reaction("COOH ethyl deprotection", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                or checker.check_reaction("TMS deprotection from alkyne", rsmi)
                or checker.check_reaction("N-glutarimide deprotection", rsmi)
                or checker.check_reaction("Phthalimide deprotection", rsmi)
            ):

                print(f"Late-stage deprotection detected at depth {depth} in reaction: {rsmi}")
                late_stage_deprotection_found = True

            # If no specific deprotection reaction was found, check for protecting groups manually
            if not late_stage_deprotection_found:
                # Check for protecting groups in reactants that are absent in product
                protecting_groups = [
                    "Boc",
                    "TMS ether protective group",
                    "Silyl protective group",
                    "Acetal/Ketal",
                ]

                for pg in protecting_groups:
                    # Check if any reactant has the protecting group
                    has_pg_in_reactants = any(
                        checker.check_fg(pg, reactant) for reactant in reactants
                    )

                    # Check if product doesn't have the protecting group
                    product_has_pg = checker.check_fg(pg, product)

                    if has_pg_in_reactants and not product_has_pg:
                        print(
                            f"Late-stage {pg} deprotection detected at depth {depth} in reaction: {rsmi}"
                        )
                        late_stage_deprotection_found = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage deprotection strategy: {late_stage_deprotection_found}")
    return late_stage_deprotection_found
