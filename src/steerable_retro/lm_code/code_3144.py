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
    This function detects if the synthetic route follows a linear synthesis strategy
    without convergent steps.

    A linear synthesis strategy has sequential reactions where each step adds to a growing
    main structure, while a convergent strategy involves separate synthesis of fragments
    that are later combined.
    """
    is_linear = True

    def is_common_reagent(smiles):
        """Identify common reagents that don't contribute significantly to the product structure"""
        # Create a list of common reagents that shouldn't count toward convergence
        common_reagents = [
            "O",
            "OO",
            "CO",
            "CCO",
            "C(=O)O",
            "C(=O)OC",
            "C(=O)OCC",
            "CN",
            "C[NH2+]",
            "C[NH3+]",
            "CC[NH2]",
            "CC[NH3+]",
            "C#N",
            "CS(=O)(=O)O",
            "CS(=O)(=O)[O-]",
            "Cl",
            "Br",
            "I",
            "F",
            "B(O)O",
            "B(OC)OC",
            "B(OCC)OCC",
            "P(OC)OC",
            "P(OCC)OCC",
            "C(=O)Cl",
            "C(=O)Br",
            "C(=O)I",
            "C(=O)F",
            "C(=O)[OH]",
            "C[N+](C)(C)C",
            "C[Si](C)(C)C",
            "C[Si](CC)(CC)CC",
            "CC(=O)O",
            "CC(=O)[O-]",
            "CC(C)(C)O",
            "CC(C)(C)OC(=O)",
            "N#N",
            "N=[N+]=[N-]",
            "[NH4+]",
            "[Na+]",
            "[K+]",
            "[Li+]",
            "S(=O)(=O)(O)O",
            "S(=O)(=O)([O-])[O-]",
        ]

        if smiles in common_reagents:
            return True

        # Check for small molecules with few atoms
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Count carbon atoms
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            # Count total heavy atoms
            heavy_atom_count = mol.GetNumHeavyAtoms()

            # Small molecules are likely reagents
            if carbon_count <= 2 and heavy_atom_count <= 8:
                return True

            # Check for common protecting groups
            if checker.check_fg("Boc", smiles) and carbon_count <= 5:
                return True
            if "Si" in smiles and carbon_count <= 4:  # Silyl protecting groups
                return True

            # Check for common boronic acids/esters
            if "B" in smiles and carbon_count <= 3:
                return True

            # Check for common coupling reagents
            if checker.check_fg("Triflate", smiles) or checker.check_fg("Tosylate", smiles):
                return True

        return False

    def get_molecular_significance(smiles):
        """Evaluate the significance of a molecule in terms of contribution to structure"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0

        score = 0

        # Carbon count contributes to significance
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        score += carbon_count

        # Rings contribute to significance
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        score += ring_count * 3  # Rings are weighted more heavily

        # Functional groups can indicate significance
        for fg in ["Ester", "Amide", "Amine", "Alcohol", "Carboxylic acid"]:
            if checker.check_fg(fg, smiles):
                score += 1

        return score

    def is_linear_reaction_type(rsmi):
        """Check if the reaction type is inherently linear"""
        linear_reaction_types = [
            "Esterification of Carboxylic Acids",
            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
            "Oxidation of aldehydes to carboxylic acids",
            "Alcohol protection with silyl ethers",
            "Alcohol deprotection from silyl ethers",
            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
            "Reduction of aldehydes and ketones to alcohols",
            "Boc amine deprotection",
            "Boc amine protection",
            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
            "Formation of NOS Heterocycles",
            "Benzimidazole formation from aldehyde",
            "Benzimidazole formation from acyl halide",
            "Benzimidazole formation from ester/carboxylic acid",
            "Benzoxazole formation from aldehyde",
            "Benzoxazole formation from acyl halide",
            "Benzoxazole formation from ester/carboxylic acid",
            "Benzoxazole formation (intramolecular)",
            "Benzothiazole formation from aldehyde",
            "Benzothiazole formation from acyl halide",
            "Benzothiazole formation from ester/carboxylic acid",
            "Paal-Knorr pyrrole synthesis",
            "Intramolecular amination (heterocycle formation)",
            "Intramolecular amination of azidobiphenyls (heterocycle formation)",
            "Pictet-Spengler",
            "Fischer indole",
            "Friedlaender chinoline",
        ]

        for reaction_type in linear_reaction_types:
            if checker.check_reaction(reaction_type, rsmi):
                print(f"Detected linear reaction type: {reaction_type}")
                return True

        # Check for cyclization reactions which are typically linear
        reactants_part = rsmi.split(">")[0]
        product_part = rsmi.split(">")[-1]

        reactant_mol = Chem.MolFromSmiles(reactants_part)
        product_mol = Chem.MolFromSmiles(product_part)

        if reactant_mol and product_mol:
            # If product has more rings than reactants, it's likely a cyclization (linear)
            reactant_rings = reactant_mol.GetRingInfo().NumRings()
            product_rings = product_mol.GetRingInfo().NumRings()

            if product_rings > reactant_rings:
                print(
                    f"Detected cyclization reaction (ring count increased from {reactant_rings} to {product_rings})"
                )
                return True

        return False

    def analyze_atom_mapping(rsmi):
        """Analyze atom mapping to determine if reaction is linear"""
        try:
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # If there's only one reactant or the reaction is a known linear type, it's likely linear
            if "." not in reactants_part or is_linear_reaction_type(rsmi):
                return True

            # Parse atom mapping from reactants and product
            reactants = reactants_part.split(".")

            # Find the reactant with the most mapped atoms that appear in the product
            max_shared_atoms = 0
            total_product_maps = 0

            # Extract atom mapping numbers from product
            product_mol = Chem.MolFromSmiles(product_part)
            if not product_mol:
                return True  # Default to linear if we can't parse the product

            product_maps = set()
            for atom in product_mol.GetAtoms():
                if atom.GetAtomMapNum() > 0:
                    product_maps.add(atom.GetAtomMapNum())

            total_product_maps = len(product_maps)
            if total_product_maps == 0:
                return True  # No atom mapping in product, default to linear

            # Check each reactant's contribution
            for reactant in reactants:
                if is_common_reagent(reactant):
                    continue

                # Extract atom mapping numbers from this reactant
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                reactant_maps = set()
                for atom in reactant_mol.GetAtoms():
                    if atom.GetAtomMapNum() > 0:
                        reactant_maps.add(atom.GetAtomMapNum())

                # Count shared mapped atoms
                shared_atoms = len(reactant_maps.intersection(product_maps))
                max_shared_atoms = max(max_shared_atoms, shared_atoms)

                # If this reactant contributes significantly to the product, it's the main reactant
                if shared_atoms > 0 and shared_atoms >= 0.5 * total_product_maps:
                    print(
                        f"Main reactant identified with {shared_atoms}/{total_product_maps} shared atoms"
                    )
                    return True

            # If no single reactant contributes significantly, it might be convergent
            if max_shared_atoms < 0.5 * total_product_maps and total_product_maps > 0:
                print(
                    f"No main reactant found, max shared atoms: {max_shared_atoms}/{total_product_maps}"
                )

                # Special case: Check if this is a heterocycle formation reaction
                # Many heterocycle formations appear convergent but are actually linear
                product_rings = product_mol.GetRingInfo().NumRings()

                # Check for common heterocycles in the product
                for ring in [
                    "benzimidazole",
                    "benzoxazole",
                    "benzothiazole",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                    "pyrrole",
                    "pyridine",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                ]:
                    if checker.check_ring(ring, product_part):
                        print(f"Detected heterocycle formation ({ring}), considering as linear")
                        return True

                return False

            return True
        except Exception as e:
            print(f"Error in atom mapping analysis: {e}")
            # Default to checking reactant count if atom mapping analysis fails
            return len(reactants_part.split(".")) <= 2

    def dfs_traverse(node):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found a convergent step
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # First check if this is a known linear reaction type
            if is_linear_reaction_type(rsmi):
                return

            # Analyze atom mapping to determine if the reaction is linear
            if not analyze_atom_mapping(rsmi):
                print(f"Atom mapping analysis suggests convergent synthesis in: {rsmi}")

                # Double-check with reactant analysis before concluding
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Filter out common reagents
                significant_reactants = []
                for reactant in reactants_smiles:
                    if not is_common_reagent(reactant):
                        significant_reactants.append(reactant)

                # If there are multiple significant reactants, check their molecular significance
                if len(significant_reactants) > 1:
                    significance_scores = [
                        get_molecular_significance(r) for r in significant_reactants
                    ]

                    # Sort reactants by significance
                    sorted_reactants = sorted(
                        zip(significant_reactants, significance_scores),
                        key=lambda x: x[1],
                        reverse=True,
                    )

                    # If the second most significant reactant has a high score relative to the most significant,
                    # it might be a convergent step
                    if len(sorted_reactants) >= 2:
                        most_significant = sorted_reactants[0][1]
                        second_significant = sorted_reactants[1][1]

                        # If the second reactant is at least 60% as significant as the first,
                        # and has a minimum significance of 6 (e.g., 2 carbons + 1 ring or 6 carbons)
                        if second_significant >= 6 and second_significant >= 0.6 * most_significant:
                            print(
                                f"Confirmed convergent step with significant reactants: {sorted_reactants}"
                            )
                            is_linear = False
                        else:
                            print(
                                f"Reactant analysis suggests this is actually linear: {sorted_reactants}"
                            )
                else:
                    print("Only one significant reactant found, considering as linear")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
