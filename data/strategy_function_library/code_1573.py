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


from rdkit import Chem
from rdkit.Chem import rdFMCS

# Refactored lists for enumeration
IMPORTANT_RINGS = [
    "benzene", "pyridine", "pyrimidine", "pyrazine", "indole", "quinoline",
    "isoquinoline", "naphthalene", "thiophene", "furan", "pyrrole", "imidazole",
    "oxazole", "thiazole", "triazole", "tetrazole", "piperidine", "piperazine",
    "morpholine", "tetrahydrofuran", "cyclohexane",
]

COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf", "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters", "Stille reaction_aryl", "Stille reaction_vinyl",
    "Stille reaction_benzyl", "Stille reaction_allyl", "Stille reaction_aryl OTf",
    "Stille reaction_vinyl OTf", "Stille reaction_benzyl OTf", "Stille reaction_allyl OTf",
    "Stille reaction_other", "Stille reaction_other OTf", "Negishi coupling",
    "Heck terminal vinyl", "Heck_terminal_vinyl", "Heck_non-terminal_vinyl",
    "Oxidative Heck reaction", "Oxidative Heck reaction with vinyl ester",
    "Heck reaction with vinyl ester and amine", "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_aryl OTf", "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl halide", "Sonogashira alkyne_alkenyl OTf",
    "Sonogashira acetylene_alkenyl OTf", "Sonogashira alkyne_acyl halide",
    "Sonogashira acetylene_acyl halide",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", "Ullmann condensation",
    "Ullmann-Goldberg Substitution amine", "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol", "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride", "Goldberg coupling", "Suzuki", "Stille",
    "decarboxylative_coupling", "Hiyama-Denmark Coupling", "Kumada cross-coupling",
    "Aryllithium cross-coupling", "Chan-Lam alcohol", "Chan-Lam amine", "Chan-Lam etherification",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) fragment coupling strategy.
    This is defined as a reaction from a specific list of coupling reactions
    (e.g., Suzuki, Stille, Buchwald-Hartwig) that joins at least two substantial
    fragments. A fragment is considered substantial if it is of a certain size
    or contains a ring system from a predefined list of important scaffolds.
    """
    late_coupling = False
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
        nonlocal late_coupling, findings_json

        # Check if this is a reaction node in the late stage (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is a known coupling reaction
                is_coupling_reaction = False
                detected_coupling_reactions = []
                for rxn_type in COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling_reaction = True
                        detected_coupling_reactions.append(rxn_type)
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

                if is_coupling_reaction:
                    # Record the positional constraint if met
                    if depth <= 1:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "coupling_reaction",
                                "position": "final_or_penultimate_step"
                            }
                        })

                    # Count substantial reactants (excluding small molecules and reagents)
                    substantial_reactants = []
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol is None:
                            continue

                        # Consider size and structural complexity
                        heavy_atoms = mol.GetNumHeavyAtoms()
                        ring_count = len(Chem.GetSSSR(mol))

                        # Check for important ring structures
                        has_important_ring = False
                        detected_rings = []
                        for ring in IMPORTANT_RINGS:
                            if checker.check_ring(ring, reactant):
                                has_important_ring = True
                                detected_rings.append(ring)
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)

                        # A substantial fragment has either:
                        # - At least 7 heavy atoms, or
                        # - At least 5 heavy atoms and contains a ring, or
                        # - Contains an important ring structure
                        if (
                            (heavy_atoms >= 7)
                            or (heavy_atoms >= 5 and ring_count > 0)
                            or has_important_ring
                        ):
                            substantial_reactants.append(reactant)

                    # If at least two substantial fragments are coupled, flag it.
                    if len(substantial_reactants) >= 2:
                        late_coupling = True
                        # Record the count constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "substantial_reactant_in_coupling_reaction",
                                "operator": ">=",
                                "value": 2
                            }
                        })

            except Exception:
                # Silently ignore errors in reaction analysis to not halt the entire process
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'chemical' or other, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)
    return late_coupling, findings_json