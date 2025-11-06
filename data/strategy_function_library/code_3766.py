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


# Refactored lists as module-level constants
COMPLEXITY_RINGS = [
    "benzene", "pyridine", "furan", "pyrrole", "thiophene", "imidazole",
    "cyclopropane", "cyclobutane", "cyclopentane", "cyclohexane", "cycloheptane",
    "pyrazole", "oxazole", "thiazole", "indole", "naphthalene",
]

COMPLEXITY_FGS = [
    "Ester", "Primary amide", "Secondary amide", "Tertiary amide",
    "Primary alcohol", "Secondary alcohol", "Tertiary alcohol",
    "Primary amine", "Secondary amine", "Tertiary amine",
    "Carboxylic acid", "Nitrile", "Aldehyde", "Ketone",
    "Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide",
]

JOINING_REACTIONS = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "Heck terminal vinyl", "Stille reaction_aryl", "Negishi coupling",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Williamson Ether Synthesis",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Esterification of Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Reductive amination with aldehyde", "Reductive amination with ketone",
    "Ugi reaction", "Michael addition", "aza-Michael addition primary",
    "aza-Michael addition secondary", "thia-Michael addition", "oxa-Michael addition",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides", "Alkylation of amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a convergent synthesis strategy where two or more complex fragments
    are joined together in a late-stage reaction.
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

    has_convergent_step = False
    max_depth = 0

    # First pass to determine the maximum depth of the synthesis tree
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)

    # Define what makes a molecule complex
    def is_complex_molecule(mol_smiles):
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check atom count (size complexity)
        atom_count = mol.GetNumAtoms()

        # Check ring count (structural complexity)
        ring_count = 0
        for ring_name in COMPLEXITY_RINGS:
            if checker.check_ring(ring_name, mol_smiles):
                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                ring_count += 1

        # Check functional group count (functional complexity)
        fg_count = 0
        for fg in COMPLEXITY_FGS:
            if checker.check_fg(fg, mol_smiles):
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                fg_count += 1

        # Define complexity thresholds
        is_complex = (atom_count >= 10 and (ring_count >= 1 or fg_count >= 2)) or atom_count >= 15

        return is_complex

    # Check if a reaction is a coupling or joining reaction
    def is_joining_reaction(rxn_smiles):
        # Common coupling/joining reactions
        for rxn_type in JOINING_REACTIONS:
            if checker.check_reaction(rxn_type, rxn_smiles):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True

        return False

    def is_convergent_step(reaction_smiles, depth):
        # Extract reactants
        try:
            reactants = reaction_smiles.split(">")[0].split(".")
            if not reactants or all(not r for r in reactants):
                return False
        except Exception:
            return False

        # Check if this is a late-stage reaction (in the first third of the synthesis)
        is_late_stage = depth <= max(1, max_depth // 3)
        if is_late_stage:
            # Add positional constraint if met
            if {"type": "positional", "details": {"target": "joining_reaction", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "joining_reaction", "position": "late_stage"}})
        else:
            return False

        # Check if there are at least two complex reactants
        complex_reactants = [r for r in reactants if is_complex_molecule(r)]

        # Check if this is a joining reaction
        is_joining = is_joining_reaction(reaction_smiles)

        is_convergent = len(complex_reactants) >= 2 and is_joining

        if is_convergent:
            # Add count constraint if met
            if {"type": "count", "details": {"target": "complex_reactants_in_joining_reaction", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactants_in_joining_reaction", "operator": ">=", "value": 2}})

        return is_convergent

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            if is_convergent_step(rsmi, depth):
                has_convergent_step = True

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return has_convergent_step, findings_json
