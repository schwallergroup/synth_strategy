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


HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyridine", "pyrrole", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole", "furan",
    "thiophene", "indole", "quinoline", "isoquinoline", "purine",
    "benzimidazole", "benzoxazole", "benzothiazole", "morpholine",
    "piperidine", "piperazine", "oxazoline", "thiazoline", "isoxazole",
    "isothiazole",
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki", "Buchwald-Hartwig", "N-arylation", "Stille", "Negishi",
    "Sonogashira", "Heck", "Ullmann", "Chan-Lam",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects convergent syntheses that assemble a final product from at least three distinct heterocyclic fragments, as defined in HETEROCYCLIC_RINGS_OF_INTEREST. The strategy is identified by the presence of a late-stage (depth <= 2) coupling reaction. A coupling is defined as either a named reaction from the COUPLING_REACTIONS_OF_INTEREST list or any reaction that joins at least two heterocyclic fragments.
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

    fragments = set()
    heterocyclic_fragments_in_coupling = set()
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal fragments, late_stage_coupling, heterocyclic_fragments_in_coupling, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                if checker.check_ring(ring, mol_smiles):
                    fragments.add(mol_smiles)
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break

        elif node["type"] == "reaction" and depth <= 2:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                is_coupling = False
                for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if not is_coupling and len(reactants) >= 2:
                    heterocyclic_reactants = 0
                    for reactant in reactants:
                        for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, reactant):
                                heterocyclic_reactants += 1
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break
                    if heterocyclic_reactants >= 2:
                        is_coupling = True

                if is_coupling:
                    late_stage_coupling = True
                    if {"type": "positional", "details": {"target": "coupling_reaction", "position": "depth <= 2"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "coupling_reaction", "position": "depth <= 2"}})
                    for reactant in reactants:
                        for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, reactant):
                                heterocyclic_fragments_in_coupling.add(reactant)
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    all_fragments = fragments.union(heterocyclic_fragments_in_coupling)

    result = len(all_fragments) >= 3 and late_stage_coupling

    if len(all_fragments) >= 3:
        if {"type": "count", "details": {"target": "heterocyclic_fragment", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "heterocyclic_fragment", "operator": ">=", "value": 3}})

    return result, findings_json