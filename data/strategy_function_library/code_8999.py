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


PYRAZOLE_FUNCTIONALIZATION_REACTIONS = [
    "Friedel-Crafts acylation",
    "Friedel-Crafts alkylation",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Suzuki coupling with boronic acids",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Aromatic iodination",
    "Aromatic bromination",
    "Aromatic chlorination",
    "Aromatic fluorination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route involves at least two modification events related to a pyrazole scaffold. A modification is defined as the formation of the pyrazole ring, the destruction of the ring, or its functionalization via a specific set of named reactions defined in PYRAZOLE_FUNCTIONALIZATION_REACTIONS (e.g., couplings, acylations, halogenations).
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

    pyrazole_modifications = 0
    pyrazole_present = False

    # Track which reactions have been processed to avoid double-counting
    processed_reactions = set()

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_modifications, pyrazole_present, findings_json

        if node["type"] == "mol":
            if checker.check_ring("pyrazole", node["smiles"]):
                pyrazole_present = True
                if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

        elif node["type"] == "reaction":
            reaction_id = node["metadata"].get("reaction_hash", "")
            if reaction_id in processed_reactions:
                return
            processed_reactions.add(reaction_id)

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                if checker.check_ring("pyrazole", product):
                    pyrazole_present = True
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

                    pyrazole_in_reactants = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )

                    if pyrazole_in_reactants:
                        # It's a functionalization if pyrazole is in reactants and product
                        functionalization_found = False
                        for rxn in PYRAZOLE_FUNCTIONALIZATION_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                pyrazole_modifications += 1
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                                functionalization_found = True
                        if functionalization_found:
                            # If a functionalization reaction was found, it implies pyrazole was present and modified
                            pass # pyrazole_modifications already incremented
                        else:
                            # If pyrazole is in reactants and product, but not a known functionalization, it's still a modification
                            # This case is implicitly handled by the general pyrazole_modifications += 1 below if no specific reaction was found.
                            pass
                    else:
                        # It's a formation if pyrazole is only in the product
                        pyrazole_modifications += 1
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                elif any(checker.check_ring("pyrazole", r) for r in reactants):
                    # It's a transformation/destruction if pyrazole is in reactants but not product
                    pyrazole_present = True
                    pyrazole_modifications += 1
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

            except (KeyError, IndexError):
                # Silently ignore errors from malformed reaction data
                pass

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = pyrazole_present and pyrazole_modifications >= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "pyrazole_modification",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
