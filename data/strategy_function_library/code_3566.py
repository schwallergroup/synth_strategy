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


DIAZINES_OF_INTEREST = ["pyrimidine", "pyrazine", "pyridazine"]
COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki", "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Stille", "Stille reaction_aryl", "Stille reaction_vinyl",
    "Negishi", "Heck", "Heck terminal vinyl", "Oxidative Heck reaction",
    "Sonogashira", "Sonogashira alkyne_aryl halide", "Sonogashira acetylene_aryl halide",
    "Buchwald-Hartwig", "N-arylation", "Ullmann", "Ullmann condensation",
    "Kumada cross-coupling", "Hiyama-Denmark Coupling", "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a pyridine-containing fragment and a diazine-containing fragment are joined via a cross-coupling reaction. The function specifically checks for diazines from the list: pyrimidine, pyrazine, pyridazine. The reaction must be one of a predefined list of named cross-coupling reactions (e.g., Suzuki, Stille, Buchwald-Hartwig). The strategy is confirmed only if two distinct reactant molecules, one containing only the pyridine fragment and the other containing only the diazine fragment, combine to form a product containing both.
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

    pyridine_molecules = set()
    diazine_molecules = set()
    coupling_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node.get("type") == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                if checker.check_ring("pyridine", smiles):
                    pyridine_molecules.add(smiles)
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                for d in DIAZINES_OF_INTEREST:
                    if checker.check_ring(d, smiles):
                        diazine_molecules.add(smiles)
                        if d not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(d)

        elif node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            is_coupling = False
            coupling_type = ""
            for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                if checker.check_reaction(rxn_type, rsmi):
                    is_coupling = True
                    coupling_type = rxn_type
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            if is_coupling:
                pyridine_reactants = []
                diazine_reactants = []
                for i, reactant in enumerate(reactants):
                    if checker.check_ring("pyridine", reactant):
                        pyridine_reactants.append(reactant)
                    for d in DIAZINES_OF_INTEREST:
                        if checker.check_ring(d, reactant):
                            diazine_reactants.append(reactant)

                product_has_pyridine = checker.check_ring("pyridine", product)
                product_has_diazine = any(checker.check_ring(d, product) for d in DIAZINES_OF_INTEREST)

                if (product_has_pyridine and product_has_diazine and
                    len(pyridine_reactants) > 0 and len(diazine_reactants) > 0):

                    pyridine_only = False
                    for reactant in pyridine_reactants:
                        if not any(checker.check_ring(d, reactant) for d in DIAZINES_OF_INTEREST):
                            pyridine_only = True
                            break

                    diazine_only = False
                    for reactant in diazine_reactants:
                        if not checker.check_ring("pyridine", reactant):
                            diazine_only = True
                            break

                    if pyridine_only and diazine_only:
                        coupling_reactions.append((rsmi, depth, coupling_type))
                        # Structural constraint: co-occurrence in product from distinct reactants
                        if {"type": "co-occurrence", "details": {"targets": ["pyridine", "any_diazine_of_interest"], "scope": "in the product of a cross-coupling reaction"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["pyridine", "any_diazine_of_interest"], "scope": "in the product of a cross-coupling reaction"}})
                        if {"type": "negation", "details": {"target": "co-occurrence of pyridine and a target diazine in a single reactant molecule"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "co-occurrence of pyridine and a target diazine in a single reactant molecule"}})

        for child in node.get("children", []):
            # New logic for depth calculation
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_detected = (
        len(pyridine_molecules) > 0 and len(diazine_molecules) > 0 and len(coupling_reactions) > 0
    )

    return strategy_detected, findings_json
