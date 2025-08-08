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
    This function detects a strategy involving late-stage nitrogen functionalization.
    """
    # Track if we found the key reactions
    n_functionalization_steps = []
    total_steps = 0

    print("Starting analysis for late-stage N-functionalization strategy")

    def dfs_traverse(node, depth=0):
        nonlocal n_functionalization_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            print(f"Processing reaction at depth {depth}")

            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")

                # Check for N-functionalization reactions
                is_n_functionalization = False

                # First verify that the product contains nitrogen
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in product_mol.GetAtoms())

                    if has_nitrogen:
                        print(f"Product contains nitrogen at depth {depth}")

                        # Check for N-alkylation reactions
                        if (
                            checker.check_reaction(
                                "N-alkylation of primary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction(
                                "N-alkylation of secondary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction("Methylation with MeI_primary", rsmi)
                            or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                            or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                            or checker.check_reaction("DMS Amine methylation", rsmi)
                            or checker.check_reaction(
                                "Eschweiler-Clarke Primary Amine Methylation", rsmi
                            )
                            or checker.check_reaction(
                                "Eschweiler-Clarke Secondary Amine Methylation", rsmi
                            )
                            or checker.check_reaction(
                                "Reductive methylation of primary amine with formaldehyde", rsmi
                            )
                            or checker.check_reaction("N-methylation", rsmi)
                            or checker.check_reaction("Alkylation of amines", rsmi)
                        ):
                            is_n_functionalization = True
                            print(f"Found N-alkylation at depth {depth}")

                            # Additional check for methylation specifically
                            if "I[CH3" in rsmi or "[CH3:10]" in rsmi:
                                print(f"Found methylation pattern in SMILES at depth {depth}")

                        # Check for N-acylation reactions
                        elif (
                            checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                rsmi,
                            )
                            or checker.check_reaction("Acylation of primary amines", rsmi)
                            or checker.check_reaction("Acylation of secondary amines", rsmi)
                            or checker.check_reaction(
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                            )
                            or checker.check_reaction(
                                "Acyl chloride with secondary amine to amide", rsmi
                            )
                            or checker.check_reaction(
                                "Carboxylic acid with primary amine to amide", rsmi
                            )
                            or checker.check_reaction("Ester with primary amine to amide", rsmi)
                            or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                            or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                            )
                        ):
                            is_n_functionalization = True
                            print(f"Found N-acylation at depth {depth}")

                        # Check for N-arylation reactions
                        elif (
                            checker.check_reaction(
                                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                            )
                            or checker.check_reaction(
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                rsmi,
                            )
                            or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                            or checker.check_reaction("{N-arylation_heterocycles}", rsmi)
                            or checker.check_reaction("Goldberg coupling", rsmi)
                            or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                        ):
                            is_n_functionalization = True
                            print(f"Found N-arylation at depth {depth}")

                        # Check for sulfonamide formation
                        elif (
                            checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                            )
                            or checker.check_reaction("{sulfon_amide}", rsmi)
                        ):
                            is_n_functionalization = True
                            print(f"Found sulfonamide formation at depth {depth}")

                        # Check for urea/thiourea formation
                        elif (
                            checker.check_reaction(
                                "Urea synthesis via isocyanate and primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Urea synthesis via isocyanate and secondary amine", rsmi
                            )
                            or checker.check_reaction("{urea}", rsmi)
                            or checker.check_reaction("{thiourea}", rsmi)
                        ):
                            is_n_functionalization = True
                            print(f"Found urea/thiourea formation at depth {depth}")

                        # Check for reductive amination
                        elif (
                            checker.check_reaction("Reductive amination with aldehyde", rsmi)
                            or checker.check_reaction("Reductive amination with ketone", rsmi)
                            or checker.check_reaction("Reductive amination with alcohol", rsmi)
                            or checker.check_reaction("{reductive amination}", rsmi)
                        ):
                            is_n_functionalization = True
                            print(f"Found reductive amination at depth {depth}")

                        # Check for amide formation
                        elif checker.check_reaction("Aminolysis of esters", rsmi):
                            is_n_functionalization = True
                            print(f"Found aminolysis of esters at depth {depth}")

                        # Manual check for N-methylation with MeI
                        elif "I[CH3" in rsmi and "[NH" in rsmi and "[N" in product_smiles:
                            is_n_functionalization = True
                            print(f"Found manual N-methylation pattern at depth {depth}")

                        # If we found a N-functionalization reaction, record its depth
                        if is_n_functionalization:
                            # Check if nitrogen is actually being functionalized by comparing reactants and products
                            reactant_has_nitrogen = False
                            for reactant_smiles in reactants_smiles:
                                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                                if reactant_mol and any(
                                    atom.GetAtomicNum() == 7 for atom in reactant_mol.GetAtoms()
                                ):
                                    reactant_has_nitrogen = True
                                    break

                            if reactant_has_nitrogen:
                                n_functionalization_steps.append(depth)
                                print(f"Confirmed N-functionalization at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Total steps: {total_steps}")
    print(f"N-functionalization steps at depths: {n_functionalization_steps}")

    # Check if N-functionalization occurs in the late stage of the synthesis
    # In retrosynthesis, lower depths correspond to later stages in forward synthesis
    if total_steps > 0 and n_functionalization_steps:
        # Consider the first half of steps as late stage (low depth in retrosynthesis)
        late_stage_threshold = total_steps // 2

        # In retrosynthesis, late stage reactions have low depth values
        late_stage_n_functionalizations = sum(
            1 for d in n_functionalization_steps if d <= late_stage_threshold
        )

        print(f"Late stage threshold: {late_stage_threshold}")
        print(f"Late stage N-functionalizations: {late_stage_n_functionalizations}")

        return late_stage_n_functionalizations > 0

    return False
