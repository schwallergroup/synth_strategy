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


NITROGEN_HETEROCYCLES_FOR_COUPLING = [
    "piperazine",
    "pyrazole",
    "imidazole",
    "triazole",
    "tetrazole",
    "pyrrole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage C-N bond formation between a halogenated benzothiazole and a specific nitrogen heterocycle. This function identifies convergent C-N coupling reactions (e.g., Buchwald-Hartwig, Ullmann) occurring in the final two steps of a synthesis. The target nitrogen heterocycles are defined in the NITROGEN_HETEROCYCLES_FOR_COUPLING list.
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

    # Track if we found the key features
    found_cn_coupling = False
    found_halogenated_benzothiazole = False
    found_nitrogen_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal found_cn_coupling, found_halogenated_benzothiazole, found_nitrogen_heterocycle, findings_json

        # Check for late-stage C-N bond formation reaction
        if node["type"] == "reaction" and depth <= 2:  # Late-stage reaction (final steps)
            # Add positional constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "C-N coupling with specific partners",
                    "position": "depth <= 2"
                }
            })

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-arylation reaction or any C-N bond formation
                is_n_arylation = False
                reaction_names = [
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "N-arylation_heterocycles",
                    "Buchwald-Hartwig"
                ]
                for r_name in reaction_names:
                    if checker.check_reaction(r_name, rsmi):
                        is_n_arylation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check reactants for halogenated benzothiazole and nitrogen heterocycle
                benzothiazole_reactant = None
                nitrogen_heterocycle_reactant = None

                for reactant in reactants:
                    # Check for benzothiazole
                    if checker.check_ring("benzothiazole", reactant):
                        if "benzothiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("benzothiazole")
                        # Check if it's halogenated
                        halogen_fgs = [
                            "Aromatic halide",
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                            "Alkenyl halide"
                        ]
                        for fg_name in halogen_fgs:
                            if checker.check_fg(fg_name, reactant):
                                found_halogenated_benzothiazole = True
                                benzothiazole_reactant = reactant
                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                break

                    # Check for any nitrogen heterocycle from the list
                    for ring in NITROGEN_HETEROCYCLES_FOR_COUPLING:
                        if checker.check_ring(ring, reactant):
                            found_nitrogen_heterocycle = True
                            nitrogen_heterocycle_reactant = reactant
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            break

                # Check if we have both required reactants
                if benzothiazole_reactant and nitrogen_heterocycle_reactant:
                    # Check if the product contains both rings
                    has_benzothiazole = checker.check_ring("benzothiazole", product)

                    has_n_heterocycle = False
                    for ring in NITROGEN_HETEROCYCLES_FOR_COUPLING:
                        if checker.check_ring(ring, product):
                            has_n_heterocycle = True
                            break

                    if has_benzothiazole and has_n_heterocycle:
                        # Check if the halogen is consumed in the reaction
                        halogen_count_reactants = sum(
                            1
                            for r in reactants
                            if (
                                checker.check_fg("Aromatic halide", r)
                                or checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Alkenyl halide", r)
                            )
                        )

                        halogen_count_product = sum(
                            1
                            for p in [product]
                            if (
                                checker.check_fg("Aromatic halide", p)
                                or checker.check_fg("Primary halide", p)
                                or checker.check_fg("Secondary halide", p)
                                or checker.check_fg("Tertiary halide", p)
                                or checker.check_fg("Alkenyl halide", p)
                            )
                        )

                        # If halogen is consumed or we detected an N-arylation reaction
                        if halogen_count_reactants > halogen_count_product or is_n_arylation:
                            found_cn_coupling = True
                            # Add co-occurrence constraint if met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "reaction_step",
                                    "targets": [
                                        "N-arylation reaction",
                                        "halogenated benzothiazole reactant",
                                        "specified nitrogen heterocycle reactant"
                                    ]
                                }
                            })

        # Continue traversal
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found all the key features
    result = found_cn_coupling and found_halogenated_benzothiazole and found_nitrogen_heterocycle

    return result, findings_json
