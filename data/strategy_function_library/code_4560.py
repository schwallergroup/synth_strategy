from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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

from pathlib import Path
root_data = Path(__file__).parent.parent

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


HETEROCYCLE_RINGS_OF_INTEREST = [
    "piperidine", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "triazole", "tetrazole",
    "indole", "quinoline", "isoquinoline", "benzimidazole", "benzoxazole",
    "benzothiazole", "piperazine", "morpholine",
]

NITROGEN_NUCLEOPHILES_OF_INTEREST = [
    "Primary amine", "Secondary amine", "Tertiary amine", "Aniline", "Azide",
    "Hydrazine", "Primary amide", "Secondary amide", "Tertiary amide", "Urea",
    "Hydrazone", "Substituted imine", "Unsubstituted imine", "Nitrile", "Isocyanate",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving at least two C-N bond formations and the presence or formation of a heterocycle, where at least two different types of nitrogen nucleophiles are used across the route. It identifies C-N bond forming reactions (e.g., N-arylation, reductive amination, amide coupling) and checks for specific heterocycles and nucleophiles defined in `HETEROCYCLE_RINGS_OF_INTEREST` and `NITROGEN_NUCLEOPHILES_OF_INTEREST`.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track key features
    cn_bond_formations = 0
    has_heterocycle = False
    nitrogen_nucleophile_types = set()  # To track different N nucleophile types

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations, has_heterocycle, nitrogen_nucleophile_types, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle presence and formation
            reactant_heterocycles = set()
            for reactant in reactants_smiles:
                for ring_name in HETEROCYCLE_RINGS_OF_INTEREST:
                    if checker.check_ring(ring_name, reactant):
                        reactant_heterocycles.add(ring_name)
                        has_heterocycle = True  # Update global flag
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)

            product_heterocycles = set()
            for ring_name in HETEROCYCLE_RINGS_OF_INTEREST:
                if checker.check_ring(ring_name, product_smiles):
                    product_heterocycles.add(ring_name)
                    has_heterocycle = True
                    if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring_name)

            # Check if new heterocycles were formed (without counting bond formations here)
            new_heterocycles = product_heterocycles - reactant_heterocycles

            # Check for nitrogen nucleophiles in reactants
            for reactant in reactants_smiles:
                for n_group in NITROGEN_NUCLEOPHILES_OF_INTEREST:
                    if checker.check_fg(n_group, reactant):
                        if n_group not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(n_group)

                        if n_group in ["Primary amine", "Secondary amine", "Tertiary amine"]:
                            nitrogen_nucleophile_types.add("aliphatic")
                        elif n_group == "Aniline":
                            nitrogen_nucleophile_types.add("aromatic")
                        elif n_group == "Azide":
                            nitrogen_nucleophile_types.add("azide")
                        elif n_group == "Hydrazine":
                            nitrogen_nucleophile_types.add("hydrazine")
                        elif n_group in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                            nitrogen_nucleophile_types.add("amide")
                        elif n_group in [
                            "Urea",
                            "Hydrazone",
                            "Substituted imine",
                            "Unsubstituted imine",
                        ]:
                            nitrogen_nucleophile_types.add("imine/urea")
                        elif n_group == "Nitrile":
                            nitrogen_nucleophile_types.add("nitrile")
                        elif n_group == "Isocyanate":
                            nitrogen_nucleophile_types.add("isocyanate")

            # Check for C-N bond formation reactions by specific reaction types
            if (
                checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                )
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                )
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                or checker.check_reaction("N-arylation_heterocycles", rsmi)
            ):
                cn_bond_formations += 1
                nitrogen_nucleophile_types.add("aromatic")
                if "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)")

            elif (
                checker.check_reaction("Alkylation of amines", rsmi)
                or checker.check_reaction("N-alkylation of primary amines with alkyl halides", rsmi)
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
            ):
                cn_bond_formations += 1
                nitrogen_nucleophile_types.add("aliphatic")
                if "Alkylation of amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alkylation of amines")

            elif (
                checker.check_reaction("reductive amination", rsmi)
                or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Reductive amination with alcohol", rsmi)
                or checker.check_reaction("{reductive amination}", rsmi)
            ):
                cn_bond_formations += 1
                nitrogen_nucleophile_types.add("aliphatic")
                if "reductive amination" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("reductive amination")

            elif (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
                or checker.check_reaction("Acyl chloride with primary amine to amide", rsmi)
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("Schotten-Baumann_amide", rsmi)
            ):
                cn_bond_formations += 1
                nitrogen_nucleophile_types.add("amide")
                if "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N")

            elif checker.check_reaction(
                "Paal-Knorr pyrrole synthesis", rsmi
            ) or checker.check_reaction("{Paal-Knorr pyrrole}", rsmi):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("aliphatic")
                if "Paal-Knorr pyrrole synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Paal-Knorr pyrrole synthesis")

            elif checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("heterocyclic")
                if "Formation of NOS Heterocycles" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Formation of NOS Heterocycles")

            elif (
                checker.check_reaction("Benzimidazole formation from aldehyde", rsmi)
                or checker.check_reaction("Benzimidazole formation from acyl halide", rsmi)
                or checker.check_reaction(
                    "Benzimidazole formation from ester/carboxylic acid", rsmi
                )
                or checker.check_reaction("{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi)
                or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
            ):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("aromatic")
                if "Benzimidazole formation from aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Benzimidazole formation from aldehyde")

            elif (
                checker.check_reaction("Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("{Huisgen_Cu-catalyzed_1,4-subst}", rsmi)
                or checker.check_reaction("{Huisgen_Ru-catalyzed_1,5_subst}", rsmi)
                or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
            ):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("azide")
                if "Huisgen alkyne-azide 1,3 dipolar cycloaddition" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Huisgen alkyne-azide 1,3 dipolar cycloaddition")

            elif (
                checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi)
                or checker.check_reaction("{tetrazole_terminal}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_1}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_2}", rsmi)
            ):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("azide")
                if "Azide-nitrile click cycloaddition to tetrazole" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Azide-nitrile click cycloaddition to tetrazole")

            elif checker.check_reaction("Pyrazole formation", rsmi) or checker.check_reaction(
                "{pyrazole}", rsmi
            ):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("hydrazine")
                if "Pyrazole formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Pyrazole formation")

            elif (
                checker.check_reaction("Urea synthesis via isocyanate and primary amine", rsmi)
                or checker.check_reaction("Urea synthesis via isocyanate and secondary amine", rsmi)
                or checker.check_reaction("{urea}", rsmi)
            ):
                cn_bond_formations += 1
                nitrogen_nucleophile_types.add("isocyanate")
                if "Urea synthesis via isocyanate and primary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Urea synthesis via isocyanate and primary amine")

            elif checker.check_reaction("{Pictet-Spengler}", rsmi):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("aliphatic")
                if "{Pictet-Spengler}" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("{Pictet-Spengler}")

            elif checker.check_reaction("{Fischer indole}", rsmi):
                cn_bond_formations += 1
                has_heterocycle = True
                nitrogen_nucleophile_types.add("hydrazine")
                if "{Fischer indole}" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("{Fischer indole}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy criteria are met
    strategy_detected = (
        cn_bond_formations >= 2 and has_heterocycle and len(nitrogen_nucleophile_types) >= 2
    )

    # Populate structural constraints if met
    if cn_bond_formations >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "cn_bond_formation",
                "operator": ">=",
                "value": 2,
                "description": "The route must contain at least two C-N bond forming reactions from a predefined list of named reactions."
            }
        })
    if has_heterocycle:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "heterocycle_presence_or_formation"
                ],
                "min_count": 1,
                "description": "The route must involve the presence or formation of at least one heterocycle from a predefined list."
            }
        })
    if len(nitrogen_nucleophile_types) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_nitrogen_nucleophile_types",
                "operator": ">=",
                "value": 2,
                "description": "The route must use at least two different types of nitrogen nucleophiles, which are categorized (e.g., 'aliphatic', 'aromatic', 'amide', 'azide')."
            }
        })

    return strategy_detected, findings_json
