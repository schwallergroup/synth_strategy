#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects if the synthesis route involves formation of a heterocycle in the final steps.
    Specifically looks for ring formation in the last 2 steps of synthesis.
    """
    late_stage_ring_formation = False

    # List of common heterocyclic rings to check
    heterocycle_rings = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "purine",
        "carbazole",
        "acridine",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of common heterocycle formation reactions
    heterocycle_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Fischer indole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "Pictet-Spengler",
        "Niementowski_quinazoline",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ring_formation

        if node["type"] == "reaction" and depth <= 2:  # Check final 3 steps to be safe
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected heterocycle formation reaction: {reaction_type}")
                        late_stage_ring_formation = True
                        return

                # Create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if all molecules were created successfully
                if product_mol and all(reactant_mols):
                    # Check for heterocycles in product that weren't in reactants
                    product_heterocycles = set()
                    reactant_heterocycles = set()

                    # Find heterocycles in product
                    for ring_name in heterocycle_rings:
                        if checker.check_ring(ring_name, product_smiles):
                            product_heterocycles.add(ring_name)

                    # Find heterocycles in reactants
                    for reactant in reactants_smiles:
                        for ring_name in heterocycle_rings:
                            if checker.check_ring(ring_name, reactant):
                                reactant_heterocycles.add(ring_name)

                    # Check if there are new heterocycles in the product
                    new_heterocycles = product_heterocycles - reactant_heterocycles
                    if new_heterocycles:
                        print(f"Detected new heterocycle(s) in product: {new_heterocycles}")
                        late_stage_ring_formation = True
                        return

                    # Count rings in reactants and product
                    reactant_ring_count = sum([len(Chem.GetSSSR(mol)) for mol in reactant_mols])
                    product_ring_count = len(Chem.GetSSSR(product_mol))

                    print(
                        f"Ring count - Reactants: {reactant_ring_count}, Product: {product_ring_count}"
                    )

                    # Check if ring count changed
                    if product_ring_count != reactant_ring_count:
                        # Check if the product contains a heterocycle
                        for ring in Chem.GetSymmSSSR(product_mol):
                            has_heteroatom = False
                            for atom_idx in ring:
                                atom = product_mol.GetAtomWithIdx(atom_idx)
                                if atom.GetSymbol() not in ["C", "H"]:
                                    has_heteroatom = True
                                    break
                            if has_heteroatom:
                                print(f"Detected heterocycle modification at depth {depth}")
                                late_stage_ring_formation = True
                                return

                    # Even if ring count didn't change, check for heterocycle modifications
                    # by examining if any heterocycle's atoms were modified
                    if product_heterocycles and product_heterocycles.intersection(
                        reactant_heterocycles
                    ):
                        # Check if the reaction involves modification of a heterocycle
                        # by looking for reactions that typically modify heterocycles
                        modification_reactions = [
                            "N-arylation",
                            "Buchwald-Hartwig",
                            "Ullmann-Goldberg Substitution",
                            "Chan-Lam",
                            "Aromatic substitution",
                            "Nucleophilic substitution",
                            "Aromatic hydroxylation",
                            "Aromatic nitration",
                            "Aromatic halogenation",
                        ]

                        for reaction_type in modification_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                # Check if the modification affects a heterocycle
                                for ring_name in product_heterocycles.intersection(
                                    reactant_heterocycles
                                ):
                                    # Get atom indices of the heterocycle in product
                                    product_indices = checker.get_ring_atom_indices(
                                        ring_name, product_smiles
                                    )
                                    if product_indices:
                                        print(
                                            f"Detected heterocycle modification reaction: {reaction_type}"
                                        )
                                        late_stage_ring_formation = True
                                        return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_ring_formation
