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
    This function detects a strategy where multiple heterocycles are formed
    sequentially in the synthesis route.
    """
    heterocycle_formations = []

    # List of heterocyclic rings to check - expanded list
    heterocycle_rings = [
        "isoxazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyrimidine",
        "pyridine",
        "furan",
        "thiophene",
        "pyrrole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "pyrrolidine",
        "azetidine",
        "aziridine",
        "oxirane",
        "thiirane",
        "dioxane",
        "dioxolane",
        "dithiolane",
        "thiazolidine",
        "oxazolidine",
        "pyrazine",
        "pyridazine",
        "carbazole",
        "acridine",
        "benzotriazole",
        "indazole",
        "pteridin",
        "dibenzofuran",
        "dibenzothiophene",
    ]

    # Expanded list of reactions that typically form heterocycles
    heterocycle_formation_rxns = [
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{1,2,4-triazole_acetohydrazide}",
        "{pyrazole}",
        "{oxadiazole}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "{Niementowski_quinazoline}",
        "{Fischer indole}",
        "{indole}",
        "{Friedlaender chinoline}",
        "{benzofuran}",
        "{benzothiophene}",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "{Paal-Knorr pyrrole}",
        "{triaryl-imidazole}",
        "{imidazole}",
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
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    ]

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        nonlocal heterocycle_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count heterocycles in reactants
                reactants_heterocycles = {}
                for reactant in reactants:
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, reactant):
                            if ring not in reactants_heterocycles:
                                reactants_heterocycles[ring] = 0
                            reactants_heterocycles[ring] += 1

                # Count heterocycles in product
                product_heterocycles = {}
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, product):
                        product_heterocycles[ring] = 1

                # Identify newly formed heterocycles
                formed_heterocycles = []
                for ring in product_heterocycles:
                    if (
                        ring not in reactants_heterocycles
                        or product_heterocycles[ring] > reactants_heterocycles[ring]
                    ):
                        formed_heterocycles.append(ring)

                # Check if this is a heterocycle formation reaction
                if formed_heterocycles:
                    # Check if it's a known heterocycle formation reaction
                    for rxn_type in heterocycle_formation_rxns:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(
                                f"Heterocycle formation reaction {rxn_type} detected at depth {depth}"
                            )
                            heterocycle_formations.append((depth, formed_heterocycles, rxn_type))
                            break
                    else:
                        # If no specific reaction matched but heterocycles were formed
                        # Check for any reaction that might form heterocycles
                        is_heterocycle_rxn = False
                        for rxn_type in [
                            "Cyclization",
                            "Ring formation",
                            "Ring closure",
                            "Condensation",
                            "Cycloaddition",
                            "Annulation",
                        ]:
                            # Since these generic types aren't in our list, we'll check for
                            # heterocycle-forming functional groups instead
                            if any(
                                checker.check_fg(fg, product)
                                for fg in ["Azide", "Nitrile", "Isocyanate", "Isothiocyanate"]
                            ):
                                is_heterocycle_rxn = True
                                break

                        if is_heterocycle_rxn or len(formed_heterocycles) > 0:
                            print(
                                f"Generic heterocycle formation detected at depth {depth}: {', '.join(formed_heterocycles)}"
                            )
                            heterocycle_formations.append((depth, formed_heterocycles, "generic"))

                # Also check for heterocycle modifications (not just formation)
                modified_heterocycles = []
                for ring in reactants_heterocycles:
                    if ring in product_heterocycles:
                        # The heterocycle is preserved but might be modified
                        # Check if functional groups were added to the heterocycle
                        for fg in [
                            "Primary amine",
                            "Secondary amine",
                            "Tertiary amine",
                            "Primary amide",
                            "Secondary amide",
                            "Tertiary amide",
                            "Nitrile",
                            "Nitro group",
                            "Halogen",
                            "Ester",
                            "Carboxylic acid",
                        ]:
                            if not any(
                                checker.check_fg(fg, reactant) for reactant in reactants
                            ) and checker.check_fg(fg, product):
                                modified_heterocycles.append(ring)
                                break

                if modified_heterocycles and not formed_heterocycles:
                    print(
                        f"Heterocycle modification detected at depth {depth}: {', '.join(modified_heterocycles)}"
                    )
                    heterocycle_formations.append((depth, modified_heterocycles, "modification"))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path + [node])

    # Start traversal
    dfs_traverse(route)

    # Check if there are at least 2 heterocycle formations
    print(f"Found {len(heterocycle_formations)} heterocycle formations")
    for depth, rings, rxn_type in heterocycle_formations:
        print(f"  Depth {depth}: Formed {', '.join(rings)} via {rxn_type}")

    # Check if the formations are sequential (within reasonable depth difference)
    if len(heterocycle_formations) >= 2:
        # Sort by depth to analyze sequence
        heterocycle_formations.sort(key=lambda x: x[0])

        # Check if any two formations are within a reasonable depth range (e.g., 5 steps)
        for i in range(len(heterocycle_formations) - 1):
            current_depth = heterocycle_formations[i][0]
            next_depth = heterocycle_formations[i + 1][0]

            # If formations are within 5 steps of each other, consider them sequential
            if abs(next_depth - current_depth) <= 5:
                print(
                    f"Sequential heterocycle formations detected at depths {current_depth} and {next_depth}"
                )
                return True

    # If we have at least 2 heterocycle formations, even if not strictly sequential
    return len(heterocycle_formations) >= 2
