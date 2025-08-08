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
    Detects if the synthesis route involves a convergent fragment coupling strategy
    where two complex fragments are joined via C-N bond formation.
    """
    convergent_coupling_found = False

    # C-N bond forming reactions to check
    cn_bond_forming_reactions = [
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Ugi reaction",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal convergent_coupling_found

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 3:  # Focus on late to mid-stage reactions
            print(f"Examining reaction at depth {depth}")
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Count number of reactant fragments
                reactants = reactants_part.split(".")
                print(f"Found {len(reactants)} reactant fragments")

                if len(reactants) >= 2:  # Multiple fragments being joined
                    try:
                        # Check if both reactants are complex
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                        product_mol = Chem.MolFromSmiles(product)

                        if all(reactant_mols) and product_mol:
                            # Check complexity based on atoms, rings, and functional groups
                            complex_fragments = []
                            for mol in reactant_mols:
                                mol_smiles = Chem.MolToSmiles(mol)
                                atom_count = mol.GetNumAtoms()
                                ring_count = len(Chem.GetSSSR(mol))

                                # Consider a fragment complex if it has >5 atoms and contains a ring
                                # or has >8 atoms regardless of rings
                                if (atom_count > 5 and ring_count > 0) or atom_count > 8:
                                    complex_fragments.append(mol)
                                    print(
                                        f"Complex fragment found: {mol_smiles}, atoms: {atom_count}, rings: {ring_count}"
                                    )

                            if len(complex_fragments) >= 2:
                                print(f"Found at least 2 complex fragments at depth {depth}")

                                # Check for C-N bond forming reactions
                                for reaction_type in cn_bond_forming_reactions:
                                    if checker.check_reaction(reaction_type, rsmi):
                                        print(f"Found C-N bond forming reaction: {reaction_type}")
                                        convergent_coupling_found = True
                                        return

                                # If no specific reaction found, check for amine and halide/carbonyl reactants
                                amine_reactants = []
                                electrophile_reactants = []

                                for i, mol in enumerate(reactant_mols):
                                    mol_smiles = Chem.MolToSmiles(mol)

                                    # Check for amine functional groups
                                    if (
                                        checker.check_fg("Primary amine", mol_smiles)
                                        or checker.check_fg("Secondary amine", mol_smiles)
                                        or checker.check_fg("Tertiary amine", mol_smiles)
                                        or checker.check_fg("Aniline", mol_smiles)
                                    ):
                                        amine_reactants.append(i)
                                        print(f"Amine found in reactant {i}: {mol_smiles}")

                                    # Check for electrophile functional groups
                                    if (
                                        checker.check_fg("Primary halide", mol_smiles)
                                        or checker.check_fg("Secondary halide", mol_smiles)
                                        or checker.check_fg("Tertiary halide", mol_smiles)
                                        or checker.check_fg("Aromatic halide", mol_smiles)
                                        or checker.check_fg("Aldehyde", mol_smiles)
                                        or checker.check_fg("Ketone", mol_smiles)
                                        or checker.check_fg("Carboxylic acid", mol_smiles)
                                        or checker.check_fg("Acyl halide", mol_smiles)
                                        or checker.check_fg("Ester", mol_smiles)
                                    ):
                                        electrophile_reactants.append(i)
                                        print(f"Electrophile found in reactant {i}: {mol_smiles}")

                                # Check if we have both amine and electrophile reactants
                                if amine_reactants and electrophile_reactants:
                                    print("Found both amine and electrophile reactants")

                                    # Check for C-N bond formation using multiple SMARTS patterns
                                    c_n_bond_patterns = [
                                        Chem.MolFromSmarts("[#6]-[#7]"),  # C-N single bond
                                        Chem.MolFromSmarts("[#6]=[#7]"),  # C=N double bond
                                        Chem.MolFromSmarts("[#6]:[#7]"),  # C:N aromatic bond
                                    ]

                                    product_cn_bonds = sum(
                                        len(product_mol.GetSubstructMatches(pattern))
                                        for pattern in c_n_bond_patterns
                                    )

                                    reactant_cn_bonds = sum(
                                        sum(
                                            len(mol.GetSubstructMatches(pattern))
                                            for pattern in c_n_bond_patterns
                                        )
                                        for mol in reactant_mols
                                    )

                                    print(
                                        f"C-N bonds in product: {product_cn_bonds}, in reactants: {reactant_cn_bonds}"
                                    )

                                    if product_cn_bonds > reactant_cn_bonds:
                                        print(f"Found new C-N bond formation at depth {depth}")
                                        convergent_coupling_found = True
                    except Exception as e:
                        print(f"Error processing SMILES in convergent coupling detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            if (
                not convergent_coupling_found
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Convergent coupling strategy found: {convergent_coupling_found}")
    return convergent_coupling_found
