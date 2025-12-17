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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


SUZUKI_REACTION_NAMES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki",
]

HETEROCYCLE_FORMATION_REACTION_NAMES = [
    "Formation of NOS Heterocycles",
    "benzimidazole_derivatives_aldehyde",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "Paal-Knorr pyrrole synthesis",
    "benzimidazole formation from aldehyde",
    "benzimidazole formation from acyl halide",
    "benzimidazole formation from ester/carboxylic acid",
    "benzoxazole formation from aldehyde",
    "benzoxazole formation from acyl halide",
    "benzoxazole formation from ester/carboxylic acid",
    "benzoxazole formation (intramolecular)",
    "benzothiazole formation from aldehyde",
    "benzothiazole formation from acyl halide",
    "benzothiazole formation from ester/carboxylic acid",
    "pyrazole formation",
    "tetrazole",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "oxadiazole",
    "imidazole",
    "Fischer indole",
    "indole",
]

HETEROCYCLES_OF_INTEREST = [
    "thiazole", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "furan", "pyran", "pyrimidine", "triazole",
    "tetrazole", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "quinoline", "isoquinoline", "purine",
    "piperidine", "piperazine", "morpholine", "thiomorpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy that combines early-stage heterocycle formation with a late-stage Suzuki coupling. It identifies these reactions by checking against curated lists of named reactions (e.g., 'Fischer indole', 'Suzuki coupling') or by structural changes, such as the formation of a new ring from the `HETEROCYCLES_OF_INTEREST` list or the co-occurrence of Suzuki-specific functional groups (boronic acid/ester and aromatic halide/triflate).
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

    has_suzuki_coupling = False
    has_heterocycle_formation = False
    suzuki_depth = -1
    heterocycle_depth = -1

    def is_suzuki_coupling(reaction_smiles):
        """Check if the reaction is a Suzuki coupling"""
        nonlocal findings_json
        for reaction_name in SUZUKI_REACTION_NAMES:
            if checker.check_reaction(reaction_name, reaction_smiles):
                if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                return True

        reactants = reaction_smiles.split(">")[0].split(".")
        has_boronic = False
        has_halide = False

        for reactant in reactants:
            if not reactant:
                continue
            if checker.check_fg("Boronic acid", reactant):
                has_boronic = True
                if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
            if checker.check_fg("Boronic ester", reactant):
                has_boronic = True
                if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")

            if checker.check_fg("Aromatic halide", reactant):
                has_halide = True
                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
            if checker.check_fg("Triflate", reactant):
                has_halide = True
                if "Triflate" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Triflate")

        if has_boronic and has_halide:
            return True

        return False

    def is_heterocycle_formation(reaction_smiles):
        """Check if the reaction forms a heterocycle"""
        nonlocal findings_json
        for reaction_name in HETEROCYCLE_FORMATION_REACTION_NAMES:
            if checker.check_reaction(reaction_name, reaction_smiles):
                if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                return True

        reactants = reaction_smiles.split(">")[0].split(".")
        products = reaction_smiles.split(">")[-1].split(".")

        for ring_name in HETEROCYCLES_OF_INTEREST:
            ring_in_reactants = False
            for reactant in reactants:
                if reactant and checker.check_ring(ring_name, reactant):
                    ring_in_reactants = True
                    break

            ring_in_products = False
            for product in products:
                if product and checker.check_ring(ring_name, product):
                    ring_in_products = True
                    if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                    break

            if ring_in_products and not ring_in_reactants:
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling, has_heterocycle_formation, suzuki_depth, heterocycle_depth, findings_json

        if node["type"] == "reaction":
            reaction_smiles = None
            if "metadata" in node:
                if "mapped_reaction_smiles" in node["metadata"]:
                    reaction_smiles = node["metadata"]["mapped_reaction_smiles"]
                elif "smiles" in node["metadata"]:
                    reaction_smiles = node["metadata"]["smiles"]

            if reaction_smiles:
                if is_suzuki_coupling(reaction_smiles):
                    has_suzuki_coupling = True
                    suzuki_depth = depth

                if is_heterocycle_formation(reaction_smiles):
                    has_heterocycle_formation = True
                    heterocycle_depth = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    strategy_present = (
        has_suzuki_coupling and has_heterocycle_formation and suzuki_depth < heterocycle_depth
    )

    if has_suzuki_coupling and has_heterocycle_formation:
        if {"type": "co-occurrence", "details": {"targets": ["suzuki_coupling", "heterocycle_formation"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["suzuki_coupling", "heterocycle_formation"]}})

    if suzuki_depth < heterocycle_depth and has_suzuki_coupling and has_heterocycle_formation:
        if {"type": "sequence", "details": {"before": "heterocycle_formation", "after": "suzuki_coupling"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "heterocycle_formation", "after": "suzuki_coupling"}})

    return strategy_present, findings_json
