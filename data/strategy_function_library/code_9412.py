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


HETEROCYCLES_OF_INTEREST = [
    "thiazole",
    "triazole",
    "pyridine",
    "pyrrolidine",
    "pyrimidine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "morpholine",
    "piperidine",
    "piperazine",
    "thiomorpholine",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "tetrazole",
    "furan",
    "pyran",
    "dioxane",
    "tetrahydrofuran",
    "tetrahydropyran",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy of linear assembly of heterocyclic fragments via
    sequential heteroatom bond formations (C-N amide, C-S thioethers)
    with late-stage O-alkylation to introduce a basic amine side chain.
    This check specifically identifies heterocycles from the HETEROCYCLES_OF_INTEREST list.
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
    thioether_formations = 0
    amide_formations = 0
    late_stage_o_alkylation_with_amine = False

    # Set to track unique heterocycles
    heterocycle_types = set()

    # Track reactions at each depth
    reactions_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal thioether_formations, amide_formations
        nonlocal late_stage_o_alkylation_with_amine, heterocycle_types, findings_json

        if node["type"] == "mol":
            # Check for heterocycles in molecule
            mol_smiles = node["smiles"]

            # Check for different heterocycle types
            for heterocycle in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(heterocycle, mol_smiles):
                    heterocycle_types.add(heterocycle)
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

        elif node["type"] == "reaction":
            if depth not in reactions_by_depth:
                reactions_by_depth[depth] = []

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Store reaction info for analysis
                reactions_by_depth[depth].append(
                    {"rsmi": rsmi, "reactants": reactants, "product": product}
                )

                # Check for thioether formation (C-S-C)
                thioether_rxn_names = [
                    "S-alkylation of thiols",
                    "S-alkylation of thiols (ethyl)",
                    "S-alkylation of thiols with alcohols",
                    "S-alkylation of thiols with alcohols (ethyl)",
                    "thioether_nucl_sub"
                ]
                thioether_rxn_detected = False
                for r_name in thioether_rxn_names:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        thioether_rxn_detected = True
                        break

                has_thiol_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Aromatic thiol", r):
                        if "Aromatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic thiol")
                        has_thiol_in_reactants = True
                    if checker.check_fg("Aliphatic thiol", r):
                        if "Aliphatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aliphatic thiol")
                        has_thiol_in_reactants = True

                has_monosulfide_in_product = checker.check_fg("Monosulfide", product)
                if has_monosulfide_in_product and "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")

                has_monosulfide_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Monosulfide", r):
                        has_monosulfide_in_reactants = True

                if thioether_rxn_detected or (
                    has_thiol_in_reactants
                    and has_monosulfide_in_product
                    and not has_monosulfide_in_reactants
                ):
                    thioether_formations += 1
                    if "thioether_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("thioether_formation")

                # Check for amide formation
                amide_rxn_names = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines"
                ]
                amide_rxn_detected = False
                for r_name in amide_rxn_names:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        amide_rxn_detected = True
                        break

                has_amine_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r):
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        has_amine_in_reactants = True
                    if checker.check_fg("Secondary amine", r):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        has_amine_in_reactants = True

                has_acyl_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Acyl halide", r):
                        if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                        has_acyl_in_reactants = True
                    if checker.check_fg("Carboxylic acid", r):
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        has_acyl_in_reactants = True
                    if checker.check_fg("Ester", r):
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        has_acyl_in_reactants = True

                has_amide_in_product = False
                if checker.check_fg("Primary amide", product):
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    has_amide_in_product = True
                if checker.check_fg("Secondary amide", product):
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    has_amide_in_product = True
                if checker.check_fg("Tertiary amide", product):
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                    has_amide_in_product = True

                has_amide_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Primary amide", r) or \
                       checker.check_fg("Secondary amide", r) or \
                       checker.check_fg("Tertiary amide", r):
                        has_amide_in_reactants = True

                if amide_rxn_detected or (
                    has_amine_in_reactants
                    and has_acyl_in_reactants
                    and has_amide_in_product
                    and not has_amide_in_reactants
                ):
                    amide_formations += 1
                    if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

                # Check for late-stage O-alkylation with amine introduction
                if depth <= 1:  # Final or penultimate step
                    o_alkylation_rxn_names = [
                        "Williamson Ether Synthesis",
                        "O-alkylation of carboxylic acids with diazo compounds",
                        "O-alkylation of amides with diazo compounds",
                        "Williamson ether"
                    ]
                    o_alkylation_rxn_detected = False
                    for r_name in o_alkylation_rxn_names:
                        if checker.check_reaction(r_name, rsmi):
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            o_alkylation_rxn_detected = True
                            break

                    has_alcohol_in_reactants = False
                    for r in reactants:
                        if checker.check_fg("Primary alcohol", r):
                            if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                            has_alcohol_in_reactants = True
                        if checker.check_fg("Secondary alcohol", r):
                            if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                            has_alcohol_in_reactants = True
                        if checker.check_fg("Tertiary alcohol", r):
                            if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                            has_alcohol_in_reactants = True
                        if checker.check_fg("Phenol", r):
                            if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                            has_alcohol_in_reactants = True

                    has_ether_in_product = checker.check_fg("Ether", product)
                    if has_ether_in_product and "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")

                    has_halide_with_amine = False
                    for r in reactants:
                        is_halide = (checker.check_fg("Primary halide", r) or
                                     checker.check_fg("Secondary halide", r) or
                                     checker.check_fg("Tertiary halide", r))
                        is_amine = (checker.check_fg("Primary amine", r) or
                                    checker.check_fg("Secondary amine", r) or
                                    checker.check_fg("Tertiary amine", r))
                        if is_halide:
                            if checker.check_fg("Primary halide", r) and "Primary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                            if checker.check_fg("Secondary halide", r) and "Secondary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                            if checker.check_fg("Tertiary halide", r) and "Tertiary halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")
                        if is_amine:
                            if checker.check_fg("Primary amine", r) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if checker.check_fg("Secondary amine", r) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            if checker.check_fg("Tertiary amine", r) and "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                        if is_halide and is_amine:
                            has_halide_with_amine = True

                    if (
                        o_alkylation_rxn_detected
                        and has_halide_with_amine
                        and has_alcohol_in_reactants
                        and has_ether_in_product
                    ):
                        late_stage_o_alkylation_with_amine = True
                        if {"type": "positional", "details": {"target": "o_alkylation_with_amine_side_chain", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "o_alkylation_with_amine_side_chain", "position": "late_stage"}})

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_present = (
        len(heterocycle_types) >= 2
        and (thioether_formations + amide_formations >= 1)
        and (late_stage_o_alkylation_with_amine or thioether_formations + amide_formations >= 2)
    )

    # Record structural constraints based on final flags
    if len(heterocycle_types) >= 2:
        if {"type": "count", "details": {"target": "unique_heterocycles_from_list", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "unique_heterocycles_from_list", "operator": ">=", "value": 2}})

    if (thioether_formations + amide_formations) >= 1:
        if {"type": "count", "details": {"target": "amide_or_thioether_formations", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "amide_or_thioether_formations", "operator": ">=", "value": 1}})

    if (thioether_formations + amide_formations) >= 2:
        if {"type": "count", "details": {"target": "amide_or_thioether_formations", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "amide_or_thioether_formations", "operator": ">=", "value": 2}})

    return strategy_present, findings_json
