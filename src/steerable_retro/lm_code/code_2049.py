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


def main(route):
    """
    Detects if the synthesis follows a linear strategy with no convergent steps.
    Each reaction should have only one reactant that contributes to the product structure.
    Reagents, catalysts, and solvents are not counted as reactants for this purpose.
    """
    is_linear = True

    def is_reagent_or_solvent(smiles):
        """
        Identifies common reagents and solvents that shouldn't count as reactants
        for determining linearity of synthesis.
        """
        # Common solvents
        solvents = [
            "O",
            "CO",
            "CCO",
            "CCCO",
            "CCCCO",
            "CC(C)O",
            "C1CCCCC1",
            "CC(=O)O",
            "CC(=O)OC",
            "CN(C)C=O",
            "CS(=O)C",
            "ClCCl",
            "ClC(Cl)Cl",
            "FC(F)F",
            "FC(F)(F)C(F)(F)F",
        ]

        # Common reagents
        reagents = [
            "O=S(Cl)Cl",
            "O=C(Cl)Cl",
            "Cl",
            "Br",
            "I",
            "F",
            "O=C=O",
            "N#N",
            "O=O",
            "[H][H]",
            "BH3",
            "B(O)O",
            "B(OC)OC",
            "B(OH)3",
            "P(Cl)Cl",
            "P(=O)(Cl)Cl",
            "S(=O)(=O)(Cl)Cl",
            "[Na+]",
            "[K+]",
            "[Li+]",
            "[Mg+2]",
            "[Zn+2]",
            "[Cu+]",
            "[Cu+2]",
            "[Pd]",
            "[Pt]",
        ]

        # Simple molecules like water, ammonia, HCl, etc.
        simple_molecules = [
            "O",
            "N",
            "Cl",
            "Br",
            "I",
            "F",
            "[H][H]",
            "CO2",
            "NH3",
            "HCl",
            "HBr",
            "HI",
            "HF",
        ]

        # Clean the SMILES by removing atom mapping
        clean_smiles = "".join([c for c in smiles if not (c.isdigit() or c == ":")])
        mol = Chem.MolFromSmiles(clean_smiles)

        if mol is None:
            return False

        # Check if it's a small molecule (likely a reagent/solvent)
        if mol.GetNumAtoms() <= 3:
            return True

        # Check against known lists
        for known_smiles in solvents + reagents + simple_molecules:
            known_mol = Chem.MolFromSmiles(known_smiles)
            if known_mol and mol.GetNumAtoms() == known_mol.GetNumAtoms():
                if mol.HasSubstructMatch(known_mol) or known_mol.HasSubstructMatch(mol):
                    return True

        return False

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Count significant reactants (excluding reagents/solvents)
            significant_reactants = [r for r in reactants if not is_reagent_or_solvent(r)]

            # If there's more than one significant reactant, it's not a linear synthesis
            if len(significant_reactants) > 1:
                is_linear = False
                print(f"Non-linear step detected: {rsmi}")
                print(f"Significant reactants: {significant_reactants}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy detected")
    else:
        print("Non-linear synthesis strategy detected")

    return is_linear
