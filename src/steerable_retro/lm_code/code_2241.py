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
    Detects linear synthesis routes that build complexity through sequential heterocycle formation.
    """
    # Track reaction steps and heterocycle formations
    reaction_count = 0
    heterocycle_formations = []  # Track formations with depth
    linear_synthesis = True

    # List of heterocyclic rings to check
    heterocyclic_rings = [
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

    # List of heterocycle formation reactions to check
    heterocycle_formation_reactions = [
        "benzothiazole_formation_from_aldehyde",
        "benzothiazole_formation_from_acyl_halide",
        "benzothiazole_formation_from_ester/carboxylic_acid",
        "benzoxazole_formation_from_aldehyde",
        "benzoxazole_formation_from_acyl_halide",
        "benzoxazole_formation_from_ester/carboxylic_acid",
        "benzoxazole_formation_(intramolecular)",
        "benzimidazole_formation_from_aldehyde",
        "benzimidazole_formation_from_acyl_halide",
        "benzimidazole_formation_from_ester/carboxylic_acid",
        "Paal-Knorr_pyrrole_synthesis",
        "Formation_of_NOS_Heterocycles",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{3-nitrile-pyridine}",
        "{pyrazole}",
        "{Paal-Knorr pyrrole}",
        "{triaryl-imidazole}",
        "{Fischer indole}",
        "{benzofuran}",
        "{benzothiophene}",
        "{indole}",
        "{oxadiazole}",
        "{imidazole}",
    ]

    # Common reagents that shouldn't count toward convergence
    common_reagents = [
        "Carbon dioxide",
        "Carbon monoxide",
        "Methanol",
        "Zinc halide",
        "Magnesium halide",
        "Alkyl lithium",
        "Aryl lithium",
    ]

    def is_reagent(smiles):
        """Check if a molecule is likely a reagent rather than a main reactant"""
        for reagent in common_reagents:
            if checker.check_fg(reagent, smiles):
                return True
        # Check for small molecules that are often reagents
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() < 3:  # Small molecules are often reagents
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_synthesis

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Filter out empty reactants
            reactants_smiles = [r for r in reactants_smiles if r.strip()]

            # Count main reactants (excluding common reagents)
            main_reactants = [r for r in reactants_smiles if not is_reagent(r)]

            # Check if this is a convergent step (more than 3 main reactants)
            if len(main_reactants) > 3:
                linear_synthesis = False
                print(f"Found convergent step with {len(main_reactants)} main reactants")

            # Check for heterocycle formation using reaction checkers
            heterocycle_formed = False
            formed_ring = None

            # First check specific reaction types
            for reaction_type in heterocycle_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    heterocycle_formed = True
                    print(f"Found heterocycle formation reaction: {reaction_type}")
                    heterocycle_formations.append((depth, reaction_type))
                    break

            # If no specific reaction type was found, check for heterocycle appearance
            if not heterocycle_formed:
                # Check if product has a heterocycle that reactants don't have
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, product_smiles):
                        # Check if any reactant has this heterocycle
                        reactants_have_ring = False
                        for reactant in reactants_smiles:
                            if checker.check_ring(ring, reactant):
                                reactants_have_ring = True
                                break

                        if not reactants_have_ring:
                            heterocycle_formed = True
                            formed_ring = ring
                            print(f"Found heterocycle formation: {ring}")
                            heterocycle_formations.append((depth, ring))
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least one heterocycle formation
    sequential_formation = len(heterocycle_formations) > 0

    # If we have multiple formations, check if they're at different stages of synthesis
    if len(heterocycle_formations) > 1:
        # Sort by depth (ascending)
        heterocycle_formations.sort(key=lambda x: x[0])
        # Check if at least some formations happen at different depths
        depths = [d for d, _ in heterocycle_formations]
        if len(set(depths)) >= 2:  # At least two different depths
            sequential_formation = True
        print(f"Heterocycle formation depths: {depths}")

    # Return True if it's a linear synthesis with at least one heterocycle formation
    # and has at least 3 reaction steps
    result = linear_synthesis and len(heterocycle_formations) > 0 and reaction_count >= 3
    print(
        f"Linear synthesis: {linear_synthesis}, Heterocycle formations: {len(heterocycle_formations)}, Reaction count: {reaction_count}, Sequential: {sequential_formation}"
    )
    return result
