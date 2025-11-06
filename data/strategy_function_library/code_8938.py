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


RING_FORMATION_REACTIONS = [
    "Diels-Alder",
    "Paal-Knorr pyrrole synthesis",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "Pictet-Spengler",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "pyrazole",
]

RING_TYPES = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "furan", "thiophene", "benzene", "naphthalene", "indole", "quinoline",
    "isoquinoline", "pyrimidine", "pyrazine", "triazole", "tetrazole",
    "cyclopropane", "cyclobutane", "cyclopentane", "cyclohexane",
    "benzimidazole", "benzoxazole", "benzothiazole", "piperidine",
    "morpholine", "pyrrolidine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses an early-stage ring formation strategy by checking for specific reactions or structural changes.

    This function identifies ring formation in any step other than the final one (i.e., depth > 1).
    A ring formation event is detected if:
    1. The reaction is a known ring-forming reaction from the RING_FORMATION_REACTIONS list.
    2. A new ring from the RING_TYPES list appears in the product.
    3. The types of rings change between reactants and products, even if the total count is the same.
    4. The total number of rings in the product is greater than in the reactants.
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

    early_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal early_ring_formation, findings_json

        if node["type"] == "reaction" and depth > 1:
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                reactant_mols = [mol for mol in reactant_mols if mol is not None]
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol or not reactant_mols:
                    return

                reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)
                product_rings = product_mol.GetRingInfo().NumRings()

                ring_formation_detected = False

                # 1. Check if this is a known ring formation reaction
                for rxn_type in RING_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        ring_formation_detected = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                # 3. Check which specific rings are formed
                formed_rings = []
                for ring_type in RING_TYPES:
                    if checker.check_ring(ring_type, product):
                        if not any(
                            checker.check_ring(ring_type, r)
                            for r in reactants
                            if Chem.MolFromSmiles(r)
                        ):
                            formed_rings.append(ring_type)
                            if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_type)

                if formed_rings:
                    ring_formation_detected = True

                # 4. Check for ring transformations (same count but different types)
                if (
                    not ring_formation_detected
                    and product_rings == reactant_rings
                    and product_rings > 0
                ):
                    product_ring_types = set()
                    reactant_ring_types = set()

                    for ring_type in RING_TYPES:
                        if checker.check_ring(ring_type, product):
                            product_ring_types.add(ring_type)
                            if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_type)

                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol and checker.check_ring(ring_type, r):
                                reactant_ring_types.add(ring_type)

                    new_ring_types = product_ring_types - reactant_ring_types
                    if new_ring_types:
                        ring_formation_detected = True

                # 5. Check for cyclization reactions that might not be in our predefined list
                if not ring_formation_detected and product_rings > reactant_rings:
                    ring_formation_detected = True
                    # Add a generic 'ring_formation' to named_reactions if not already present
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                if ring_formation_detected:
                    early_ring_formation = True
                    # Add the structural constraint if early ring formation is detected
                    # This constraint is defined in the input JSON as:
                    # {"type": "positional", "details": {"target": "ring_formation", "position": "not_last_stage"}}
                    constraint_obj = {"type": "positional", "details": {"target": "ring_formation", "position": "not_last_stage"}}
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)

            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return early_ring_formation, findings_json