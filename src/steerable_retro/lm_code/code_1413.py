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
    This function detects orthogonal protection-deprotection strategy with Boc and Cbz groups
    on a piperazine scaffold or similar nitrogen-containing heterocycles.
    """
    # Track protection/deprotection events
    boc_protection = False
    boc_deprotection = False
    cbz_protection = False
    cbz_deprotection = False
    n_heterocycle_present = False

    def dfs_traverse(node):
        nonlocal boc_protection, boc_deprotection, cbz_protection, cbz_deprotection, n_heterocycle_present

        if node["type"] == "mol":
            # Check if molecule contains piperazine or similar N-heterocycles
            if (
                checker.check_ring("piperazine", node["smiles"])
                or checker.check_ring("morpholine", node["smiles"])
                or checker.check_ring("diazepane", node["smiles"])
                or checker.check_ring("piperidine", node["smiles"])
                or checker.check_ring("pyrrolidine", node["smiles"])
                or checker.check_ring("azepane", node["smiles"])
                or checker.check_ring("azetidine", node["smiles"])
                or checker.check_ring("aziridine", node["smiles"])
            ):
                n_heterocycle_present = True
                print(f"Detected N-heterocycle in: {node['smiles']}")

        elif node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for Boc protection
                if not checker.check_fg("Boc", reactants_smiles) and checker.check_fg(
                    "Boc", product_smiles
                ):
                    # Check for any Boc protection reaction
                    if any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Boc amine protection",
                            "Boc amine protection explicit",
                            "Boc amine protection with Boc anhydride",
                            "Boc amine protection (ethyl Boc)",
                            "Boc amine protection of secondary amine",
                            "Boc amine protection of primary amine",
                        ]
                    ):
                        boc_protection = True
                        print(f"Detected Boc protection: {rsmi}")
                    # Fallback check if reaction type not detected
                    elif "C(C)(C)OC(=O)" in product_smiles:
                        boc_protection = True
                        print(f"Detected Boc protection (pattern match): {rsmi}")

                # Check for Boc deprotection
                if checker.check_fg("Boc", reactants_smiles) and not checker.check_fg(
                    "Boc", product_smiles
                ):
                    if any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Boc amine deprotection",
                            "Boc amine deprotection of guanidine",
                            "Boc amine deprotection to NH-NH2",
                            "Tert-butyl deprotection of amine",
                        ]
                    ):
                        boc_deprotection = True
                        print(f"Detected Boc deprotection: {rsmi}")
                    # Fallback check if reaction type not detected
                    elif "C(C)(C)OC(=O)" in reactants_smiles:
                        boc_deprotection = True
                        print(f"Detected Boc deprotection (pattern match): {rsmi}")

                # Check for Cbz protection
                if not checker.check_fg("Carbamic ester", reactants_smiles) and checker.check_fg(
                    "Carbamic ester", product_smiles
                ):
                    # Check for benzyl pattern in the product
                    if "c1ccccc1" in product_smiles and "OCc1ccccc1" in product_smiles:
                        cbz_protection = True
                        print(f"Detected Cbz protection: {rsmi}")

                # Check for Cbz deprotection
                if checker.check_fg("Carbamic ester", reactants_smiles) and not checker.check_fg(
                    "Carbamic ester", product_smiles
                ):
                    # Check for benzyl pattern in the reactants
                    if "c1ccccc1" in reactants_smiles and "OCc1ccccc1" in reactants_smiles:
                        if any(
                            checker.check_reaction(rxn, rsmi)
                            for rxn in [
                                "Carboxyl benzyl deprotection",
                                "Hydrogenolysis of amides/imides/carbamates",
                            ]
                        ):
                            cbz_deprotection = True
                            print(f"Detected Cbz deprotection: {rsmi}")
                        # Fallback check if reaction type not detected
                        else:
                            cbz_deprotection = True
                            print(f"Detected Cbz deprotection (pattern match): {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have N-heterocycle and evidence of orthogonal protection
    # This means we need both Boc and Cbz groups involved in protection/deprotection
    strategy_present = n_heterocycle_present and (
        (boc_protection or boc_deprotection) and (cbz_protection or cbz_deprotection)
    )

    print(f"N-heterocycle present: {n_heterocycle_present}")
    print(f"Boc protection: {boc_protection}")
    print(f"Boc deprotection: {boc_deprotection}")
    print(f"Cbz protection: {cbz_protection}")
    print(f"Cbz deprotection: {cbz_deprotection}")

    if strategy_present:
        print("Orthogonal protection-deprotection strategy detected")
    else:
        print("No orthogonal protection-deprotection strategy detected")

    return strategy_present
