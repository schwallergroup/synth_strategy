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
    This function detects a strategy involving a late-stage cyclization
    (final step involves ring formation).
    """
    # Track if the final step involves ring formation
    final_step_is_cyclization = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_cyclization

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction step at depth {depth}: {rsmi}")

                # Check if metadata explicitly marks this as a ring-forming reaction
                if node.get("metadata", {}).get("RingBreaker", False):
                    print(f"Detected RingBreaker flag in metadata")
                    final_step_is_cyclization = True
                    return

                # Check if this is a known cyclization reaction
                cyclization_reaction_types = [
                    "Formation of NOS Heterocycles",
                    "Intramolecular transesterification/Lactone formation",
                    "Paal-Knorr pyrrole synthesis",
                    "Intramolecular amination (heterocycle formation)",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                    "Benzothiazole formation from aldehyde",
                    "Benzothiazole formation from acyl halide",
                    "Benzothiazole formation from ester/carboxylic acid",
                    "Benzoxazole formation from aldehyde",
                    "Benzoxazole formation from acyl halide",
                    "Benzoxazole formation from ester/carboxylic acid",
                    "Benzoxazole formation (intramolecular)",
                    "Benzimidazole formation from aldehyde",
                    "Benzimidazole formation from acyl halide",
                    "Benzimidazole formation from ester/carboxylic acid",
                    "Mitsunobu aryl ether (intramolecular)",
                    "Pictet-Spengler",
                    "Fischer indole",
                    "Diels-Alder",
                    "Diels-Alder (ON bond)",
                    "Pauson-Khand reaction",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                    "Pyrazole formation",
                    "Michael-induced ring closure from hydrazone",
                    "Michael-induced ring closure from diazoalkane",
                    "[3+2]-cycloaddition of hydrazone and alkyne",
                    "[3+2]-cycloaddition of hydrazone and alkene",
                    "[3+2]-cycloaddition of diazoalkane and alkyne",
                    "[3+2]-cycloaddition of diazoalkane and alkene",
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
                    "Azide-nitrile click cycloaddition to tetrazole",
                    "Azide-nitrile click cycloaddition to triazole",
                    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                ]

                for rxn_type in cyclization_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected cyclization reaction: {rxn_type}")
                        final_step_is_cyclization = True
                        return

                # Convert to RDKit molecules for further analysis
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(reactants):
                    # Check for ring formation by comparing ring counts
                    product_rings = Chem.GetSSSR(product)
                    reactant_rings_total = sum(Chem.GetSSSR(r) for r in reactants)

                    if len(product_rings) > reactant_rings_total:
                        print(
                            f"Detected ring formation: Product has {len(product_rings)} rings, reactants have {reactant_rings_total} rings"
                        )
                        final_step_is_cyclization = True
                        return

                    # Check for specific ring types in the product that might not be in reactants
                    ring_types = [
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
                        "trioxane",
                        "dioxepane",
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
                        "aziridine",
                        "azetidine",
                        "azepane",
                        "diazepane",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "purine",
                        "carbazole",
                        "acridine",
                        "thiophene",
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
                        "benzoxazole",
                        "benzothiazole",
                        "benzimidazole",
                        "pteridin",
                        "phenothiazine",
                        "phenoxazine",
                        "dibenzofuran",
                        "dibenzothiophene",
                        "xanthene",
                        "thioxanthene",
                        "pyrroline",
                        "pyrrolidone",
                        "imidazolidine",
                        "porphyrin",
                        "indazole",
                        "benzotriazole",
                    ]

                    for ring_type in ring_types:
                        # Check if ring exists in product but not in all reactants
                        if checker.check_ring(ring_type, product_smiles):
                            print(f"Found {ring_type} in product")
                            ring_in_all_reactants = all(
                                checker.check_ring(ring_type, r) for r in reactants_smiles
                            )
                            if not ring_in_all_reactants:
                                print(f"Detected {ring_type} formation in step at depth {depth}")
                                final_step_is_cyclization = True
                                return

                    # Check for intramolecular reactions (often cyclizations)
                    if len(reactants_smiles) == 1 and "." not in reactants_smiles[0]:
                        # One reactant becoming one product often indicates intramolecular reaction
                        reactant_mol = reactants[0]

                        # Check if any functional groups that often participate in cyclization are present
                        cyclization_prone_fgs = [
                            "Carboxylic acid",
                            "Ester",
                            "Amide",
                            "Amine",
                            "Alcohol",
                            "Alkyne",
                            "Alkene",
                            "Azide",
                            "Nitrile",
                            "Isocyanate",
                        ]

                        fg_count_reactant = sum(
                            1
                            for fg in cyclization_prone_fgs
                            if checker.check_fg(fg, reactants_smiles[0])
                        )
                        fg_count_product = sum(
                            1
                            for fg in cyclization_prone_fgs
                            if checker.check_fg(fg, product_smiles)
                        )

                        if fg_count_reactant > fg_count_product:
                            print(
                                f"Detected potential intramolecular cyclization: functional group count decreased"
                            )
                            final_step_is_cyclization = True
                            return

            except Exception as e:
                print(f"Error analyzing reaction step: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final result: late_stage_cyclization = {final_step_is_cyclization}")
    return final_step_is_cyclization
