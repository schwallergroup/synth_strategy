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
    This function detects a synthetic strategy involving N-demethylation.
    """
    n_demethylation_detected = False

    def dfs_traverse(node):
        nonlocal n_demethylation_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # First check if this is explicitly an N-demethylation reaction
                if checker.check_reaction("N-demethylation", rsmi):
                    print(f"N-demethylation reaction detected directly: {rsmi}")
                    n_demethylation_detected = True
                    return

                # Skip if this is a methylation reaction (opposite of demethylation)
                methylation_reactions = [
                    "N-methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "Parnes methylation",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "DMS Amine methylation",
                ]

                for rxn_type in methylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Methylation reaction detected ({rxn_type}), not demethylation")
                        return

                # In retrosynthesis, the product is our starting material and reactants are targets
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if not all(reactant_mols) or not product_mol:
                    print("Warning: Could not parse some molecules")
                    return

                # Check for nitrogen-containing functional groups in reactants and product
                n_containing_fgs = [
                    "Primary amine",
                    "Secondary amine",
                    "Tertiary amine",
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                    "Aniline",
                    "Pyridine",
                    "Pyrrole",
                    "Imidazole",
                    "Triazole",
                ]

                for i, reactant_smi in enumerate(reactants_smiles):
                    if not reactant_smi:
                        continue

                    # Check for N-containing functional groups in reactant
                    reactant_n_fgs = []
                    for fg in n_containing_fgs:
                        if checker.check_fg(fg, reactant_smi):
                            reactant_n_fgs.append(fg)

                    if not reactant_n_fgs:
                        continue

                    # Check for N-containing functional groups in product
                    product_n_fgs = []
                    for fg in n_containing_fgs:
                        if checker.check_fg(fg, product_smiles):
                            product_n_fgs.append(fg)

                    if not product_n_fgs:
                        continue

                    # Check for demethylation pattern
                    # In retrosynthesis: secondary→tertiary or primary→secondary indicates demethylation strategy
                    demethylation_patterns = [
                        ("Secondary amine", "Tertiary amine"),
                        ("Primary amine", "Secondary amine"),
                        ("Secondary amide", "Tertiary amide"),
                        ("Primary amide", "Secondary amide"),
                    ]

                    for reactant_fg, product_fg in demethylation_patterns:
                        if reactant_fg in reactant_n_fgs and product_fg in product_n_fgs:
                            print(
                                f"Potential N-demethylation pattern detected: {reactant_fg} → {product_fg}"
                            )

                            # Additional check: verify by comparing molecular formulas
                            reactant_mol = reactant_mols[i]

                            # Calculate formula differences
                            reactant_formula = AllChem.CalcMolFormula(reactant_mol)
                            product_formula = AllChem.CalcMolFormula(product_mol)

                            print(
                                f"Reactant formula: {reactant_formula}, Product formula: {product_formula}"
                            )

                            # Check if product has more carbon and hydrogen atoms (consistent with methyl group)
                            # This is a simplification - a more robust approach would track atom mappings
                            reactant_c_count = sum(
                                1 for a in reactant_mol.GetAtoms() if a.GetSymbol() == "C"
                            )
                            product_c_count = sum(
                                1 for a in product_mol.GetAtoms() if a.GetSymbol() == "C"
                            )

                            if product_c_count > reactant_c_count:
                                # In retrosynthesis, product (starting material) should have more carbon atoms
                                # than reactant (target) if a methyl group was removed
                                print(
                                    f"N-demethylation confirmed: product C count {product_c_count} > reactant C count {reactant_c_count}"
                                )
                                n_demethylation_detected = True
                                return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return n_demethylation_detected
