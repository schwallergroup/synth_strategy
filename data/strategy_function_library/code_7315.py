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


from rdkit import Chem
from rdkit.Chem import Descriptors

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Heck terminal vinyl",
    "Oxidative Heck reaction",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille reaction_aryl OTf",
    "Stille reaction_vinyl OTf",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage convergent coupling strategy, defined as a final reaction step (depth=0) that joins at least two complex fragments. A fragment's complexity is determined by its size and ring count. The reaction must be one of the types specified in the COUPLING_REACTIONS_OF_INTEREST list.
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

    convergent_final_step = False

    def is_complex_fragment(smiles):
        """Check if a molecule is a complex fragment based on size and structure"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mw = Descriptors.MolWt(mol)
        ring_count = mol.GetNumRings()
        heavy_atom_count = mol.GetNumHeavyAtoms()

        complexity_score = 0
        if mw > 150:
            complexity_score += 1
        if ring_count >= 1:
            complexity_score += 1
        if heavy_atom_count > 10:
            complexity_score += 1

        return complexity_score >= 2

    def dfs_traverse(node, depth=0):
        nonlocal convergent_final_step, findings_json

        if node["type"] == "mol" and depth == 0:
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    dfs_traverse(child, depth)
            return

        if node["type"] == "reaction" and depth == 0:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                if len(reactants_smiles) >= 2:
                    complex_reactants = [s for s in reactants_smiles if is_complex_fragment(s)]

                    is_coupling = False
                    for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

                    if is_coupling and len(complex_reactants) >= 2:
                        convergent_final_step = True
                        # Add structural constraints if conditions are met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "coupling_reaction",
                                "position": "last_stage"
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants",
                                "scope": "last_stage_reaction",
                                "operator": ">=",
                                "value": 2
                            }
                        })
            except (KeyError, IndexError):
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # i.e., if current node is 'mol'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return convergent_final_step, findings_json
