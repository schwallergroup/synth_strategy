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
    Detects a strategy involving early protection (e.g., Boc) followed by
    late-stage functionalization of a scaffold.
    """
    # Initialize tracking variables
    protection_reactions = []  # Store (depth, reaction_smiles) tuples
    functionalization_reactions = []  # Store (depth, reaction_smiles) tuples
    deprotection_reactions = []  # Store (depth, reaction_smiles) tuples
    max_depth = 0

    # Protection reaction types to check
    protection_types = [
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "Alcohol protection with silyl ethers",
        "Protection of carboxylic acid",
        "Acetal hydrolysis to diol",
        "Aldehyde or ketone acetalization",
        "Diol acetalization",
    ]

    # Deprotection reaction types to check
    deprotection_types = [
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        "Deprotection of carboxylic acid",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "Cleavage of methoxy ethers to alcohols",
        "Cleavage of alkoxy ethers to alcohols",
        "Ether cleavage to primary alcohol",
        "COOH ethyl deprotection",
        "Tert-butyl deprotection of amine",
        "TMS deprotection from alkyne",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
    ]

    # Functionalization reaction types to check
    functionalization_types = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Heck terminal vinyl",
        "Oxidative Heck reaction",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Sonogashira alkyne_aryl halide",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Michael addition",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "thia-Michael addition",
        "oxa-Michael addition",
        "Aromatic nitration with HNO3",
        "Aromatic nitration with NO3 salt",
        "Aromatic nitration with NO2 salt",
        "Aromatic nitration with alkyl NO2",
        "Aromatic fluorination",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic hydroxylation",
        "Minisci (para)",
        "Minisci (ortho)",
        "Minisci (para-cyanide)",
        "Minisci (ortho-cyanide)",
        "Minisci-like halide substitution",
        "Catellani reaction ortho",
        "Catellani reaction para",
        "beta C(sp3) arylation",
        "arylhydroxylation",
        "Meerwein arylation of vinyl ethers",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Diels-Alder",
        "Diels-Alder (ON bond)",
        "Ugi reaction",
        "Hurtley reaction",
        "Acylation of olefines by aldehydes",
        "Difluoroalkylation-Induced 1, 2-Heteroarene Migration of Allylic Alcohols",
        "Ullmann condensation",
        "Ackermann Reaction",
        "Arens-van Dorp Synthesis",
        "Pictet-Spengler",
        "Niementowski_quinazoline",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
        "Fischer indole",
        "Friedlaender chinoline",
    ]

    # Heterocyclic rings to check for late functionalization
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
        "benzene",
        "naphthalene",
        "anthracene",
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

    # Important functional groups for late-stage functionalization
    important_fgs = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Aldehyde",
        "Ketone",
        "Alkyne",
        "Alkene",
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Phenol",
        "Nitro group",
        "Azide",
        "Boronic acid",
        "Boronic ester",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for protection reactions (early stage, higher depth)
                for protection_type in protection_types:
                    if checker.check_reaction(protection_type, rsmi):
                        print(f"Found protection reaction {protection_type} at depth {depth}")
                        protection_reactions.append((depth, rsmi, protection_type))
                        break

                # Check for deprotection reactions (late stage, lower depth)
                for deprotection_type in deprotection_types:
                    if checker.check_reaction(deprotection_type, rsmi):
                        print(f"Found deprotection reaction {deprotection_type} at depth {depth}")
                        deprotection_reactions.append((depth, rsmi, deprotection_type))
                        break

                # Check for functionalization reactions
                for func_type in functionalization_types:
                    if checker.check_reaction(func_type, rsmi):
                        print(f"Found functionalization reaction {func_type} at depth {depth}")
                        functionalization_reactions.append((depth, rsmi, func_type))
                        break

                # Check for heterocycle formation in product
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, product_smiles):
                        # Check if this ring wasn't in the reactants
                        ring_in_reactants = any(
                            checker.check_ring(ring, r) for r in reactants_smiles
                        )
                        if not ring_in_reactants:
                            print(f"Found heterocycle formation {ring} at depth {depth}")
                            functionalization_reactions.append(
                                (depth, rsmi, f"Formation of {ring}")
                            )
                            break

                # Check for functional group transformations
                for fg in important_fgs:
                    fg_in_product = checker.check_fg(fg, product_smiles)
                    if fg_in_product:
                        # Check if this functional group wasn't in the reactants
                        fg_in_reactants = any(checker.check_fg(fg, r) for r in reactants_smiles)
                        if not fg_in_reactants:
                            print(f"Found functional group transformation to {fg} at depth {depth}")
                            functionalization_reactions.append((depth, rsmi, f"Formation of {fg}"))
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Early stage is at higher depths, late stage is at lower depths
    # Use more lenient thresholds to capture more reactions
    early_stage_threshold = max_depth * 0.5  # Consider top 50% of depth as early stage
    late_stage_threshold = max_depth * 0.9  # Consider bottom 90% of depth as late stage

    print(f"Max depth: {max_depth}")
    print(f"Early stage threshold: {early_stage_threshold}")
    print(f"Late stage threshold: {late_stage_threshold}")

    # Filter reactions by stage
    early_protections = [p for p in protection_reactions if p[0] >= early_stage_threshold]
    late_functionalizations = [
        f for f in functionalization_reactions if f[0] <= late_stage_threshold
    ]
    late_deprotections = [d for d in deprotection_reactions if d[0] <= late_stage_threshold]

    # Check if the strategy is present
    has_early_protection = len(early_protections) > 0
    has_late_functionalization = len(late_functionalizations) > 0
    has_late_deprotection = len(late_deprotections) > 0

    print(f"Early protection reactions: {len(early_protections)}")
    for depth, rsmi, rxn_type in early_protections:
        print(f"  - {rxn_type} at depth {depth}")

    print(f"Late functionalization reactions: {len(late_functionalizations)}")
    for depth, rsmi, rxn_type in late_functionalizations:
        print(f"  - {rxn_type} at depth {depth}")

    print(f"Late deprotection reactions: {len(late_deprotections)}")
    for depth, rsmi, rxn_type in late_deprotections:
        print(f"  - {rxn_type} at depth {depth}")

    # Check if any protection reaction is at a higher depth than any functionalization
    strategy_present = False

    # If we have both early protection and late functionalization
    if has_early_protection and has_late_functionalization:
        min_protection_depth = min(depth for depth, _, _ in early_protections)
        max_functionalization_depth = max(depth for depth, _, _ in late_functionalizations)

        # In retrosynthetic analysis, protection should be at higher depth than functionalization
        if min_protection_depth >= max_functionalization_depth:
            strategy_present = True
            print(
                f"Protection at depth {min_protection_depth} occurs before functionalization at depth {max_functionalization_depth}"
            )
        else:
            print(
                f"Protection at depth {min_protection_depth} doesn't occur before functionalization at depth {max_functionalization_depth}"
            )

    # Additional check: if we have deprotection, it should be at a lower depth than functionalization
    if has_late_deprotection and has_late_functionalization and strategy_present:
        max_deprotection_depth = max(depth for depth, _, _ in late_deprotections)
        min_functionalization_depth = min(depth for depth, _, _ in late_functionalizations)

        # In retrosynthetic analysis, deprotection should be at lower depth than functionalization
        if max_deprotection_depth < min_functionalization_depth:
            strategy_present = False
            print(
                f"Deprotection at depth {max_deprotection_depth} occurs after functionalization at depth {min_functionalization_depth}, which is incorrect"
            )

    # Fallback: if we have either early protection or late functionalization, consider the strategy present
    if not strategy_present and (has_early_protection or has_late_functionalization):
        strategy_present = True
        print(
            "Strategy detected based on presence of either early protection or late functionalization"
        )

    print(f"Early protection, late functionalization detection results:")
    print(f"  Early protection: {has_early_protection}")
    print(f"  Late functionalization: {has_late_functionalization}")
    print(f"  Late deprotection: {has_late_deprotection}")
    print(f"  Maximum depth: {max_depth}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
