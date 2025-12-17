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
from rdkit.Chem import rdMolDescriptors

CONVERGENT_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Buchwald-Hartwig",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Heck terminal vinyl",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling",
    "Wittig reaction with triphenylphosphorane",
    "Wittig with Phosphonium",
    "Grignard from aldehyde to alcohol",
    "Grignard from ketone to alcohol",
    "Mitsunobu aryl ether",
    "Diels-Alder",
    "Michael addition",
    "Aldol condensation",
    "Ugi reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy. A reaction is flagged as convergent if it occurs in the latter half of the synthesis (depth <= max_depth / 2), joins at least two 'complex' fragments (defined as having > 5 heavy atoms or at least one ring), and is one of the specified reaction types in the CONVERGENT_REACTION_TYPES list (e.g., Suzuki, Buchwald-Hartwig, Wittig, Grignard, etc.).
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

    is_convergent = False
    max_depth = 0

    # First pass to determine the maximum depth of the synthesis
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)

    # Second pass to detect convergent synthesis
    def dfs_traverse(node, depth, max_depth):
        nonlocal is_convergent, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a reaction with multiple reactants in the second half of synthesis
            if len(reactants) > 1 and all(r.strip() for r in reactants) and depth <= (max_depth // 2 + 1):
                is_fragment_joining_reaction = False
                for rxn_type in CONVERGENT_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_fragment_joining_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_fragment_joining_reaction:
                    # Record positional constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "convergent_reaction",
                            "position": "latter_half"
                        }
                    })

                    # Count complex reactants (significant fragments)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            num_heavy_atoms = mol.GetNumHeavyAtoms()
                            num_rings = rdMolDescriptors.CalcNumRings(mol)

                            # Consider reactants with significant structure as complex
                            if num_heavy_atoms >= 6 or num_rings > 0:
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        is_convergent = True
                        # Record count constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants_in_convergent_reaction",
                                "operator": ">=",
                                "value": 2
                            }
                        })

        # Continue traversing
        for child in node.get("children", []):
            if is_convergent:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases for reaction child
                new_depth = depth + 1
            # If current node is reaction, depth remains the same for chemical child
            
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal from the root node at depth 0. The first reaction will be at depth 1.
    dfs_traverse(route, 0, max_depth)
    return is_convergent, findings_json
