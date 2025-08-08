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
    Detects if the synthesis uses a late-stage Suzuki coupling strategy where two complex fragments
    are connected in the final step via a Suzuki coupling, with one fragment being prepared through
    borylation and the other containing a carboxylic acid derivative.
    """
    # Initialize tracking variables
    final_suzuki_coupling = False
    suzuki_coupling_depth = float("inf")
    borylation_present = False
    carboxylic_acid_derivative = False

    def dfs_traverse(node, depth=0):
        nonlocal final_suzuki_coupling, suzuki_coupling_depth, borylation_present, carboxylic_acid_derivative

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Suzuki coupling (prioritize depth 0 or 1 as final step)
            # Check for Suzuki coupling using the checker function
            suzuki_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "{Suzuki}",
                "Suzuki",
            ]

            for suzuki_type in suzuki_types:
                if checker.check_reaction(suzuki_type, rsmi):
                    final_suzuki_coupling = True
                    if depth < suzuki_coupling_depth:
                        suzuki_coupling_depth = depth
                    print(f"Detected {suzuki_type} at depth {depth}")
                    break

            # If not detected by reaction checker, try to check for components
            if not final_suzuki_coupling:
                has_boronic_compound = False
                has_halide_or_triflate = False

                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic_compound = True
                        print(f"Found boronic compound in reaction: {reactant}")

                    if (
                        checker.check_fg("Aromatic halide", reactant)
                        or checker.check_fg("Triflate", reactant)
                        or checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                        or checker.check_fg("Tertiary halide", reactant)
                        or "Br" in reactant
                        or "I" in reactant
                        or "Cl" in reactant
                    ):
                        has_halide_or_triflate = True
                        print(f"Found halide/triflate in reaction: {reactant}")

                # Check if product has a new biaryl bond
                if has_boronic_compound and has_halide_or_triflate:
                    final_suzuki_coupling = True
                    if depth < suzuki_coupling_depth:
                        suzuki_coupling_depth = depth
                    print(f"Detected Suzuki coupling pattern at depth {depth}")

                    # Try to identify if a new biaryl bond was formed
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            for bond in prod_mol.GetBonds():
                                a1 = bond.GetBeginAtom()
                                a2 = bond.GetEndAtom()
                                if (
                                    a1.GetIsAromatic()
                                    and a2.GetIsAromatic()
                                    and a1.GetSymbol() == "C"
                                    and a2.GetSymbol() == "C"
                                ):
                                    print("Confirmed biaryl bond formation in product")
                                    break
                    except Exception as e:
                        print(f"Error checking biaryl formation: {e}")

            # Check for borylation (typically at depth 1 or greater)
            if depth >= 1:
                # Check for borylation reactions using the checker function
                borylation_types = [
                    "Preparation of boronic acids",
                    "Preparation of boronic acids without boronic ether",
                    "Preparation of boronic acids from trifluoroborates",
                    "Preparation of boronic ethers",
                    "Synthesis of boronic acids",
                ]

                for borylation_type in borylation_types:
                    if checker.check_reaction(borylation_type, rsmi):
                        borylation_present = True
                        print(f"Detected {borylation_type} at depth {depth}")
                        break

                # If not detected by reaction checker, check for conversion of halide to boronic compound
                if not borylation_present:
                    has_halide_reactant = False
                    has_boronic_product = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or "Br" in reactant
                            or "I" in reactant
                            or "Cl" in reactant
                        ):
                            has_halide_reactant = True
                            print(f"Found halide reactant in borylation: {reactant}")

                    if checker.check_fg("Boronic acid", product) or checker.check_fg(
                        "Boronic ester", product
                    ):
                        has_boronic_product = True
                        print(f"Found boronic product in borylation: {product}")

                    if has_halide_reactant and has_boronic_product:
                        borylation_present = True
                        print(f"Detected borylation pattern at depth {depth}")

            # Check for carboxylic acid derivative transformations
            carboxylic_derivatives = [
                "Carboxylic acid",
                "Ester",
                "Acyl halide",
                "Anhydride",
                "Primary amide",
                "Secondary amide",
                "Tertiary amide",
            ]

            for fg in carboxylic_derivatives:
                for reactant in reactants:
                    if checker.check_fg(fg, reactant):
                        carboxylic_acid_derivative = True
                        print(f"Detected {fg} in reactant: {reactant}")

                if checker.check_fg(fg, product):
                    carboxylic_acid_derivative = True
                    print(f"Detected {fg} in product: {product}")

            # Also check for reactions that produce carboxylic acid derivatives
            carboxylic_reactions = [
                "Oxidation of aldehydes to carboxylic acids",
                "Oxidation of alcohol to carboxylic acid",
                "Oxidation of ketone to carboxylic acid",
                "Oxidation of nitrile to carboxylic acid",
                "Aryl halide to carboxylic acid",
                "Oxidation of amide to carboxylic acid",
                "Esterification of Carboxylic Acids",
                "Carboxylic acid to amide conversion",
            ]

            for rxn in carboxylic_reactions:
                if checker.check_reaction(rxn, rsmi):
                    carboxylic_acid_derivative = True
                    print(f"Detected {rxn} reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if all criteria are met
    # For late-stage Suzuki, it should be at depth 0 or 1
    strategy_present = (
        final_suzuki_coupling
        and borylation_present
        and carboxylic_acid_derivative
        and suzuki_coupling_depth <= 1
    )

    print(f"Late-stage Suzuki coupling strategy detected: {strategy_present}")
    print(
        f"- Final Suzuki coupling: {final_suzuki_coupling} at depth {suzuki_coupling_depth if final_suzuki_coupling else 'N/A'}"
    )
    print(f"- Borylation present: {borylation_present}")
    print(f"- Carboxylic acid derivative: {carboxylic_acid_derivative}")

    return strategy_present
