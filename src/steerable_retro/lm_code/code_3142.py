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
    This function detects if the synthetic route involves formation of a heterocycle
    in the final step (late-stage heterocycle formation).
    """
    final_heterocycle_formation = False

    # List of common heterocycles to check
    heterocycle_types = [
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
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
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
    ]

    # List of heterocycle-forming reactions
    heterocycle_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzothiazole formation from aldehyde",
        "benzothiazole formation from acyl halide",
        "benzothiazole formation from ester/carboxylic acid",
        "benzoxazole formation from aldehyde",
        "benzoxazole formation from acyl halide",
        "benzoxazole formation from ester/carboxylic acid",
        "benzoxazole formation (intramolecular)",
        "benzimidazole formation from aldehyde",
        "benzimidazole formation from acyl halide",
        "benzimidazole formation from ester/carboxylic acid",
        "tetrazole formation",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "{Niementowski_quinazoline}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{Huisgen_Ru-catalyzed_1,5_subst}",
        "{Huisgen_disubst-alkyne}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{3-nitrile-pyridine}",
        "{pyrazole}",
        "{Paal-Knorr pyrrole}",
        "{triaryl-imidazole}",
        "{Fischer indole}",
        "{Friedlaender chinoline}",
        "{benzofuran}",
        "{benzothiophene}",
        "{indole}",
        "{oxadiazole}",
        "{imidazole}",
        "1,2,4-oxadiazol-5(2H)-one synthesis from nitrile, hydrogen carbonate, and hydroxylamine",
    ]

    def find_final_reaction(node, depth=0, path=None):
        """Find the final reaction in the synthetic route (closest to target)"""
        if path is None:
            path = []

        # If this is the target molecule (root of the tree)
        if depth == 0 and node["type"] == "mol":
            # Look for its parent reaction
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    return child

        # If this is a reaction node
        if node["type"] == "reaction":
            # Check if this reaction produces a molecule that has no further reactions
            # (i.e., it's the final step before the target)
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    # Check if this molecule has no reaction children
                    has_reaction_children = False
                    for grandchild in child.get("children", []):
                        if grandchild["type"] == "reaction":
                            has_reaction_children = True
                            break

                    if not has_reaction_children:
                        return node

        # If not found at this level, check children
        for child in node.get("children", []):
            result = find_final_reaction(child, depth + 1, path + [node])
            if result:
                return result

        return None

    # Find the final reaction in the route
    final_reaction = find_final_reaction(route)

    if final_reaction and "metadata" in final_reaction and "rsmi" in final_reaction["metadata"]:
        print(f"Found final reaction: {final_reaction['metadata']['rsmi']}")

        rsmi = final_reaction["metadata"]["rsmi"]

        try:
            # Extract reactants and product
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants_smiles = parts[0].split(".")
                product_smiles = parts[2]

                print(f"Final reaction - Reactants: {reactants_smiles}, Product: {product_smiles}")

                # Convert to RDKit molecules
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol and all(mol is not None for mol in reactants_mols):
                    # Check if the reaction is a known heterocycle-forming reaction
                    for rxn_type in heterocycle_forming_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected heterocycle formation reaction: {rxn_type}")
                            final_heterocycle_formation = True
                            break

                    if not final_heterocycle_formation:
                        # Check for heterocycle presence in product but not in reactants
                        heterocycles_in_product = []
                        for heterocycle in heterocycle_types:
                            if checker.check_ring(heterocycle, product_smiles):
                                heterocycles_in_product.append(heterocycle)
                                print(f"Found {heterocycle} in product")

                        # For each heterocycle in the product, check if it's new
                        for heterocycle in heterocycles_in_product:
                            # Check if this heterocycle was not present in any of the reactants
                            heterocycle_in_reactants = any(
                                checker.check_ring(heterocycle, r) for r in reactants_smiles
                            )

                            if not heterocycle_in_reactants:
                                print(f"Detected new {heterocycle} formation in final step")
                                final_heterocycle_formation = True
                                break

                    # If still not identified, check for ring transformations
                    if not final_heterocycle_formation:
                        # Get all rings in reactants and product
                        reactant_rings = []
                        for r_mol in reactants_mols:
                            if r_mol:
                                ring_info = r_mol.GetRingInfo()
                                for idx in range(ring_info.NumRings()):
                                    ring_atoms = ring_info.AtomRings()[idx]
                                    reactant_rings.append(set(ring_atoms))

                        product_rings = []
                        if product_mol:
                            ring_info = product_mol.GetRingInfo()
                            for idx in range(ring_info.NumRings()):
                                ring_atoms = ring_info.AtomRings()[idx]
                                product_rings.append(set(ring_atoms))

                        # Count rings in reactants and product
                        reactants_ring_count = sum(
                            [mol.GetRingInfo().NumRings() for mol in reactants_mols if mol]
                        )
                        product_ring_count = (
                            product_mol.GetRingInfo().NumRings() if product_mol else 0
                        )

                        print(
                            f"Ring count - Reactants: {reactants_ring_count}, Product: {product_ring_count}"
                        )

                        # Check if product has different ring structure
                        if product_ring_count != reactants_ring_count:
                            # Check if product contains any heterocycle
                            for heterocycle in heterocycle_types:
                                if checker.check_ring(heterocycle, product_smiles):
                                    print(
                                        f"Detected {heterocycle} in product with changed ring count"
                                    )
                                    final_heterocycle_formation = True
                                    break

                        # Even if ring count is the same, check for ring transformations
                        # by looking for heterocycles in the product that weren't in reactants
                        if not final_heterocycle_formation:
                            # Check for specific ring-forming reactions
                            ring_forming_reactions = [
                                "Intramolecular amination",
                                "Cyclization",
                                "Ring formation",
                                "Ring closure",
                            ]

                            for reaction_pattern in ring_forming_reactions:
                                if reaction_pattern.lower() in rsmi.lower():
                                    # Check if product contains any heterocycle
                                    for heterocycle in heterocycle_types:
                                        if checker.check_ring(heterocycle, product_smiles):
                                            print(
                                                f"Detected {heterocycle} in product via {reaction_pattern}"
                                            )
                                            final_heterocycle_formation = True
                                            break
                                    if final_heterocycle_formation:
                                        break
        except Exception as e:
            print(f"Error processing reaction SMILES: {e}")
    else:
        print("No final reaction found in the route")

    print(f"Final result: {final_heterocycle_formation}")
    return final_heterocycle_formation
