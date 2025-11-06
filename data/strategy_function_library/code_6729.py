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


# Refactoring for Enumeration: Isolate lists of chemical entities
SUGAR_RING_TYPES = ['pyran', 'furan', 'tetrahydropyran', 'tetrahydrofuran']
ALCOHOL_FG_TYPES = ['Primary alcohol', 'Secondary alcohol', 'Tertiary alcohol', 'Phenol', 'Aromatic alcohol']

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a glycosylation strategy using trichloroacetimidate activation.
    It looks for:
    1. Presence of trichloroacetimidate group in early steps
    2. Formation of glycosidic bond (any glycoside, not just aryl)
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

    has_trichloroacetimidate = False
    forms_glycosidic_bond = False
    trichloroacetimidate_reactions = []
    glycosidic_bond_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal has_trichloroacetimidate, forms_glycosidic_bond, findings_json

        # In retrosynthesis, early steps have higher depth values
        is_early_step = depth >= 2

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for trichloroacetimidate group in reactants
                for reactant in reactants:
                    try:
                        # Check for trichloro group and sugar structure
                        if is_early_step and checker.check_fg("Trichloro group", reactant):
                            findings_json["atomic_checks"]["functional_groups"].append("Trichloro group")
                            # Check for sugar rings
                            has_sugar = False
                            for ring in SUGAR_RING_TYPES:
                                if checker.check_ring(ring, reactant):
                                    has_sugar = True
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)

                            # Check for C=N-O-C pattern characteristic of trichloroacetimidate
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and has_sugar:
                                pattern = Chem.MolFromSmarts("[CX3](=[NX2])[OX2][C]")
                                if mol.HasSubstructMatch(pattern):
                                    has_trichloroacetimidate = True
                                    trichloroacetimidate_reactions.append(rsmi)
                                    print(
                                        f"Found trichloroacetimidate group in reactant at depth {depth}: {reactant}"
                                    )
                                    # Record positional constraint if found in early step
                                    if {"type": "positional", "details": {"target": "trichloroacetimidate_presence", "position": "early_stage"}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "trichloroacetimidate_presence", "position": "early_stage"}})

                    except Exception as e:
                        print(f"Error checking for trichloroacetimidate: {e}")
                        continue

                # Check for glycosidic bond formation
                try:
                    # Check for sugar-like rings in product
                    has_sugar_ring = False
                    for ring in SUGAR_RING_TYPES:
                        if checker.check_ring(ring, product):
                            has_sugar_ring = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    # Check for ether formation
                    has_ether = checker.check_fg("Ether", product)
                    if has_ether and "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")

                    # Check if any reactant has alcohol (potential glycosyl acceptor)
                    has_alcohol_reactant = False
                    for r in reactants:
                        for fg in ALCOHOL_FG_TYPES:
                            if checker.check_fg(fg, r):
                                has_alcohol_reactant = True
                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)

                    # Check if trichloroacetimidate is consumed (not in product)
                    trichloroacetimidate_consumed = any(
                        checker.check_fg("Trichloro group", r) for r in reactants
                    ) and not checker.check_fg("Trichloro group", product)

                    # Glycosidic bond formation: sugar ring + ether + alcohol reactant + trichloroacetimidate consumed
                    if (
                        has_sugar_ring
                        and has_ether
                        and has_alcohol_reactant
                        and trichloroacetimidate_consumed
                    ):
                        forms_glycosidic_bond = True
                        glycosidic_bond_reactions.append(rsmi)
                        if "glycosidic_bond_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("glycosidic_bond_formation")
                        print(f"Detected glycosidic bond formation at depth {depth}")
                except Exception as e:
                    print(f"Error checking for glycosidic bond: {e}")

        # Process children nodes with increased depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    print(f"Has trichloroacetimidate: {has_trichloroacetimidate}")
    print(f"Forms glycosidic bond: {forms_glycosidic_bond}")

    result = has_trichloroacetimidate and forms_glycosidic_bond

    if result:
        print(f"Trichloroacetimidate reactions: {trichloroacetimidate_reactions}")
        print(f"Glycosidic bond reactions: {glycosidic_bond_reactions}")
        # Record co-occurrence constraint if both conditions are met
        if {"type": "co-occurrence", "details": {"targets": ["has_trichloroacetimidate", "forms_glycosidic_bond"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["has_trichloroacetimidate", "forms_glycosidic_bond"]}})

    return result, findings_json