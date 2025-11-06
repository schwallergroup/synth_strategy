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


COUPLING_REACTIONS = [
    "Suzuki",
    "Stille",
    "Stille reaction_aryl",
    "Negishi",
    "Heck",
    "Heck_terminal_vinyl",
    "Sonogashira",
    "Sonogashira alkyne_aryl halide",
    "Buchwald-Hartwig",
    "N-arylation",
    "Ullmann",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies a late-stage convergent synthesis strategy, defined as a reaction at depth 0 or 1 where at least two complex fragments are joined via a key coupling reaction. A fragment is considered 'complex' if it has at least 8 heavy atoms or a complexity score (heavy atoms + 2 * rings) of 10 or more. The reaction is identified as a key coupling if it matches a named reaction from the `COUPLING_REACTIONS` list, or if it involves the coupling of a boronic acid/ester with a halide.
    """
    is_convergent = False
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json

        if node["type"] == "reaction" and depth <= 1:
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            complex_fragments = []
            for reactant in reactants:
                if reactant.strip():
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            heavy_atom_count = mol.GetNumHeavyAtoms()
                            ring_count = mol.GetRingInfo().NumRings()
                            complexity_score = heavy_atom_count + (ring_count * 2)
                            if heavy_atom_count >= 8 or complexity_score >= 10:
                                complex_fragments.append(reactant)
                    except Exception:
                        continue

            if len(complex_fragments) >= 2:
                # Structural constraint: count complex_reactants
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "complex_reactants",
                        "operator": ">=",
                        "value": 2
                    }
                })

                is_coupling = False

                for reaction_name in COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if not is_coupling:
                    has_boronic = False
                    has_halide = False
                    for reactant in reactants:
                        if checker.check_fg("Boronic acid", reactant):
                            has_boronic = True
                            findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                        if checker.check_fg("Boronic ester", reactant):
                            has_boronic = True
                            findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")

                        if checker.check_fg("Aromatic halide", reactant):
                            has_halide = True
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        if checker.check_fg("Primary halide", reactant):
                            has_halide = True
                            findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                        if checker.check_fg("Secondary halide", reactant):
                            has_halide = True
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                        if checker.check_fg("Tertiary halide", reactant):
                            has_halide = True
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")

                    if has_boronic and has_halide:
                        is_coupling = True
                        # Structural constraint: co-occurrence boronic_acid_or_ester and halide
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "boronic_acid_or_ester",
                                    "halide"
                                ],
                                "scope": "reaction_reactants"
                            }
                        })

                if is_coupling:
                    is_convergent = True
                    # Structural constraint: positional key_coupling_reaction
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "key_coupling_reaction",
                            "position": "depth <= 1"
                        }
                    })

        for child in node.get("children", []):
            if not is_convergent:
                new_depth = depth
                if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., it's a chemical)
                    new_depth += 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    # Ensure unique entries in lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))
    # Structural constraints are objects, so we need a custom way to deduplicate if necessary
    # For this problem, we assume appending the full object is sufficient and duplicates are fine if they represent multiple instances of the same constraint being met.

    return is_convergent, findings_json