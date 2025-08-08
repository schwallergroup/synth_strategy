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
    This function detects if the synthetic route follows a linear strategy (no convergent steps).

    A linear strategy means each reaction step has only one significant reactant that contributes
    to the growing main molecule, with other reactants being reagents or small building blocks.
    """
    is_linear = True

    # List of reaction types that are considered linear despite having multiple significant reactants
    linear_compatible_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with sulfonic esters",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_alkenyl OTf",
        "Sonogashira acetylene_acyl halide",
        "Sonogashira alkyne_acyl halide",
        "Heck terminal vinyl",
        "Heck non-terminal vinyl",
        "Oxidative Heck reaction",
        "Oxidative Heck reaction with vinyl ester",
        "Heck reaction with vinyl ester and amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Negishi coupling",
        "Stille reaction_vinyl",
        "Stille reaction_aryl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_vinyl OTf",
        "Stille reaction_aryl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Aromatic bromination",
        "Aromatic chlorination",
        "Aromatic fluorination",
        "Aromatic iodination",
        "Wohl-Ziegler bromination benzyl primary",
        "Wohl-Ziegler bromination benzyl secondary",
        "Wohl-Ziegler bromination benzyl tertiary",
        "Wohl-Ziegler bromination allyl primary",
        "Wohl-Ziegler bromination allyl secondary",
        "Wohl-Ziegler bromination allyl tertiary",
    ]

    # List of functional groups that indicate significant building blocks even in small molecules
    significant_fgs = [
        "Boronic acid",
        "Boronic ester",
        "Magnesium halide",
        "Zinc halide",
        "Tin",
        "Triflate",
        "Aromatic halide",
    ]

    # Common reagents that should not be counted as significant reactants
    common_reagents = [
        "O=C1OCCCC1",  # THF
        "CN(C)C=O",  # DMF
        "CC(=O)OC",  # Methyl acetate
        "CC(C)=O",  # Acetone
        "CCO",  # Ethanol
        "CO",  # Methanol
        "CCCCN(C)C",  # DIPEA
        "c1ccncc1",  # Pyridine
        "CN1CCCC1=O",  # NMP
        "CS(=O)(=O)C",  # DMSO
        "ClCCl",  # DCM
        "ClC(Cl)Cl",  # Chloroform
        "CC(C)(C)O",  # tert-butanol
        "c1ccccc1",  # Benzene
        "Cc1ccccc1",  # Toluene
        "O",  # Water
        "N",  # Ammonia
        "CC#N",  # Acetonitrile
        "[Na+]",  # Sodium ion
        "[K+]",  # Potassium ion
        "[Li+]",  # Lithium ion
        "[Mg+2]",  # Magnesium ion
        "[Zn+2]",  # Zinc ion
        "[Cu+]",  # Copper ion
        "[Pd]",  # Palladium
        "[Pt]",  # Platinum
        "[Cl-]",  # Chloride ion
        "[Br-]",  # Bromide ion
        "[I-]",  # Iodide ion
        "[OH-]",  # Hydroxide ion
        "O=C=O",  # Carbon dioxide
        "CC1(C)C(=O)N(Br)C(=O)N1Br",  # NBS (N-Bromosuccinimide)
        "CC1(C)C(=O)N(Cl)C(=O)N1Cl",  # NCS (N-Chlorosuccinimide)
        "BrBr",  # Bromine
        "ClCl",  # Chlorine
        "II",  # Iodine
        "FC(F)(F)C(=O)OC(=O)C(F)(F)F",  # TFAA
        "BrCCBr",  # 1,2-dibromoethane
        "ClCCCl",  # 1,2-dichloroethane
    ]

    def is_reagent(smiles):
        """Check if a molecule is likely a reagent rather than a significant reactant"""
        try:
            # Check against common reagent list
            for reagent in common_reagents:
                if reagent == smiles:
                    return True

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Small molecules are often reagents
            if mol.GetNumHeavyAtoms() <= 4:
                return True

            # Check for common solvent/reagent patterns
            if (
                smiles.count("C") + smiles.count("c") + smiles.count("O") + smiles.count("N")
                == len(smiles)
                and mol.GetNumHeavyAtoms() < 8
            ):
                return True

            # Check for brominating agents
            if "N(Br)" in smiles or "NBr" in smiles:
                if "C(=O)N" in smiles and "C(=O)N" in smiles:
                    print(f"Identified brominating agent: {smiles}")
                    return True

            # Check for halogenating agents in general
            if mol.GetNumHeavyAtoms() < 15:
                if smiles.count("Br") >= 1 and "C(=O)" in smiles:
                    print(f"Possible brominating agent: {smiles}")
                    return True

            return False
        except:
            return False

    def dfs_traverse(node):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found a convergent step
            return

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a known linear-compatible reaction type
            is_linear_compatible = False
            for rxn_type in linear_compatible_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found linear-compatible reaction: {rxn_type}")
                    is_linear_compatible = True
                    break

            # Check specifically for bromination reactions
            if not is_linear_compatible:
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                # Check if the reaction is a bromination
                if product_mol and "Br" in product:
                    reactants = rsmi.split(">")[0].split(".")
                    for reactant in reactants:
                        if "Br" in reactant and (
                            "N(Br)" in reactant or "NBr" in reactant or "BrBr" in reactant
                        ):
                            print(
                                f"Identified bromination reaction with brominating agent: {reactant}"
                            )
                            is_linear_compatible = True
                            break

            if not is_linear_compatible:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count significant reactants (non-reagents)
                significant_reactants = 0
                for reactant in reactants:
                    try:
                        if is_reagent(reactant):
                            print(f"Identified reagent: {reactant}")
                            continue

                        mol = Chem.MolFromSmiles(reactant)
                        # Consider a reactant significant if it's large or contains important functional groups
                        if mol:
                            has_significant_fg = any(
                                checker.check_fg(fg, reactant) for fg in significant_fgs
                            )
                            is_large = mol.GetNumHeavyAtoms() > 8

                            if is_large or has_significant_fg:
                                significant_reactants += 1
                                print(
                                    f"Significant reactant found: {reactant} (large: {is_large}, has significant FG: {has_significant_fg})"
                                )
                    except Exception as e:
                        print(f"Error processing reactant {reactant}: {e}")
                        continue

                # If more than one significant reactant, it's a convergent step
                if significant_reactants > 1:
                    print(
                        f"Detected convergent step with {significant_reactants} significant reactants: {rsmi}"
                    )
                    is_linear = False
                else:
                    print(f"Linear step with {significant_reactants} significant reactant(s)")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return is_linear
