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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a morpholine group is introduced in a late-stage step (depth <= 2), either by direct incorporation of a morpholine-containing reactant via specific C-N bond-forming reactions (e.g., N-alkylation, N-arylation, sulfonamide formation) or by de novo ring formation.
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

    morpholine_added = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_added, late_stage, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late stage (expanded definition)
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains morpholine
                if checker.check_ring("morpholine", product):
                    findings_json["atomic_checks"]["ring_systems"].append("morpholine")

                    # Check if morpholine is being incorporated (present in reactants as separate entity)
                    morpholine_in_reactants = False
                    non_morpholine_in_reactants = False

                    for reactant in reactants:
                        if checker.check_ring("morpholine", reactant):
                            # If the reactant is primarily just morpholine (or a simple derivative)
                            if len(Chem.MolFromSmiles(reactant).GetAtoms()) <= 15:
                                morpholine_in_reactants = True
                        else:
                            non_morpholine_in_reactants = True

                    # Check for relevant reaction types that would incorporate morpholine
                    is_incorporation_reaction = False

                    # Check common reaction types for morpholine incorporation
                    if checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    ):
                        is_incorporation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")
                    elif checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    ):
                        is_incorporation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of primary amines with alkyl halides")
                    elif checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    ):
                        is_incorporation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of secondary amines with alkyl halides")
                    elif checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    ):
                        is_incorporation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append("Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine")
                    elif checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    ):
                        is_incorporation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append("Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine")

                    # If we have morpholine in reactants, non-morpholine in reactants, and it's an incorporation reaction
                    if morpholine_in_reactants and non_morpholine_in_reactants and is_incorporation_reaction:
                        morpholine_added = True
                        late_stage = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "morpholine_introduction_reaction", "position": "late_stage (depth <= 2)"}})

                    # Alternative detection: if product has morpholine but no reactant has it
                    elif not morpholine_in_reactants and non_morpholine_in_reactants:
                        # This means morpholine was formed in the reaction
                        morpholine_added = True
                        late_stage = True
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "morpholine_introduction_reaction", "position": "late_stage (depth <= 2)"}})

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return morpholine_added and late_stage, findings_json