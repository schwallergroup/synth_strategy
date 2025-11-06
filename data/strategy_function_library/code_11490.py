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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving the formation of an oxazolidinone ring
    from a protected amino alcohol.

    In the forward direction, this involves cyclization of a protected amino alcohol
    to form an oxazolidinone. In retrosynthesis (how we traverse the tree), we look
    for reactions where an oxazolidinone is converted to an amino alcohol.
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

    has_oxazolidinone_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_oxazolidinone_formation, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_has_oxazolidinone = False
                for reactant in reactants:
                    if checker.check_ring("oxazolidine", reactant) and checker.check_fg(
                        "Carbamic ester", reactant
                    ):
                        reactant_has_oxazolidinone = True
                        if "oxazolidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxazolidine")
                        if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        break

                if reactant_has_oxazolidinone:
                    product_has_oxazolidinone = checker.check_ring(
                        "oxazolidine", product
                    ) and checker.check_fg("Carbamic ester", product)

                    if not product_has_oxazolidinone:
                        has_amino = checker.check_fg("Primary amine", product) or checker.check_fg(
                            "Secondary amine", product
                        )
                        if checker.check_fg("Primary amine", product) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", product) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                        has_protected_amino = (
                            checker.check_fg("Boc", product)
                            or checker.check_fg("Carbamic ester", product)
                            or checker.check_fg("Carbamic acid", product)
                        )
                        if checker.check_fg("Boc", product) and "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                        if checker.check_fg("Carbamic ester", product) and "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        if checker.check_fg("Carbamic acid", product) and "Carbamic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic acid")

                        has_alcohol = checker.check_fg(
                            "Primary alcohol", product
                        ) or checker.check_fg("Secondary alcohol", product)
                        if checker.check_fg("Primary alcohol", product) and "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                        if checker.check_fg("Secondary alcohol", product) and "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")

                        if (has_amino or has_protected_amino) and has_alcohol:
                            has_oxazolidinone_formation = True
                            if {"type": "co-occurrence", "details": {"targets": ["ring_formation", "oxazolidine", "Carbamic ester", "Primary amine", "Secondary amine", "Boc", "Carbamic acid", "Primary alcohol", "Secondary alcohol"]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "oxazolidine", "Carbamic ester", "Primary amine", "Secondary amine", "Boc", "Carbamic acid", "Primary alcohol", "Secondary alcohol"]}})
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception as e:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return has_oxazolidinone_formation, findings_json
