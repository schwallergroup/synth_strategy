from typing import Tuple, Dict, List
import copy
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


# Refactored lists as module-level constants
HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyridine", "pyrrole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "furan", "thiophene", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole", "quinoline", "isoquinoline", "purine",
]

CN_BOND_FORMING_REACTIONS = [
    "N-arylation", "Buchwald-Hartwig", "Ullmann-Goldberg Substitution amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves a late-stage (depth <= 2) attachment of a heterocycle to a scaffold.
    This is identified by checking for a specific C-N bond forming reaction (from the CN_BOND_FORMING_REACTIONS list)
    where one reactant contains a heterocycle (from the HETEROCYCLIC_RINGS_OF_INTEREST list) and at least one other reactant does not.
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

    late_stage_hetero_attachment = False

    def is_heterocycle(mol_smiles):
        """Check if molecule contains one of the specified heterocyclic rings."""
        found_rings = []
        for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
            if checker.check_ring(ring, mol_smiles):
                found_rings.append(ring)
        if found_rings:
            # Add found rings to findings_json, ensuring no duplicates
            for r in found_rings:
                if r not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(r)
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_hetero_attachment, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late-stage reaction
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_cn_reaction = False
                found_cn_reactions = []
                for rxn in CN_BOND_FORMING_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_cn_reaction = True
                        found_cn_reactions.append(rxn)
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                if is_cn_reaction:
                    # Add positional constraint if depth <= 2
                    if depth <= 2:
                        positional_constraint = {
                            "type": "positional",
                            "details": {
                                "target": "heterocycle_attachment_reaction",
                                "condition": "depth <= 2"
                            }
                        }
                        if positional_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint)

                    # Identify heterocycle and non-heterocycle reactants
                    heterocycle_reactants = [r for r in reactants if is_heterocycle(r)]
                    non_heterocycle_reactants = [r for r in reactants if not is_heterocycle(r)]

                    # Check if one reactant is a heterocycle and another is not
                    if heterocycle_reactants and non_heterocycle_reactants:
                        # Verify the product has the heterocycle attached
                        if is_heterocycle(product):
                            late_stage_hetero_attachment = True
                            # Add co-occurrence constraint
                            co_occurrence_constraint = {
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "reaction_step",
                                    "targets": [
                                        "CN_bond_formation_reaction",
                                        "reactant_with_specified_heterocycle",
                                        "reactant_without_specified_heterocycle"
                                    ]
                                }
                            }
                            if co_occurrence_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(co_occurrence_constraint)

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            new_depth = depth + 1

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_hetero_attachment, findings_json
