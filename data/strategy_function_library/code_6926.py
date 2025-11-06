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

HETEROCYCLE_RINGS_OF_INTEREST = [
    "furan",
    "pyrrole",
    "thiophene",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "indole",
    "benzofuran",
    "benzothiophene",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
    "triazole",
    "tetrazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "purine",
    "carbazole",
    "acridine",
    "morpholine",
    "piperidine",
    "piperazine",
    "pyrrolidine",
    "tetrahydrofuran",
    "tetrahydropyran",
    "dioxane",
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
    "Stille reaction_aryl",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann condensation",
    "Heck terminal vinyl",
    "Negishi coupling",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Goldberg coupling",
    "Stille reaction_vinyl",
    "Catellani reaction ortho",
    "Catellani reaction para",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
    "Suzuki coupling with sulfonic esters",
    "Chan-Lam alcohol",
    "Chan-Lam amine",
    "Chan-Lam etherification",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy involving the coupling of two complex heterocyclic fragments
    in a convergent manner at a late stage.
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

    found_convergent_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_coupling, findings_json

        if node["type"] == "reaction" and depth <= 3:  # Late stage reaction (expanded depth)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if it's a coupling reaction from the predefined list
                is_coupling = False
                for rxn in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn, rsmi):
                        is_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                
                if is_coupling and len(reactants) >= 2:
                    # Add structural constraint for late_stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "convergent_coupling_reaction",
                            "position": "late_stage (depth <= 3)"
                        }
                    })

                    # Check if at least two reactants have heterocycles and are complex
                    heterocycle_counts = []
                    complex_fragments = []

                    for reactant in reactants:
                        if not reactant:
                            continue

                        # Count heterocycles in this reactant
                        count = 0
                        for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, reactant):
                                count += 1
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                        heterocycle_counts.append(count)

                        # Check if it's a complex fragment (>= 8 atoms)
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() >= 8:
                            complex_fragments.append(True)
                        else:
                            complex_fragments.append(False)

                    # Need at least 2 reactants with heterocycles and complexity
                    num_heterocyclic_reactants = sum(1 for count in heterocycle_counts if count >= 1)
                    num_complex_reactants = sum(complex_fragments)

                    if (len(heterocycle_counts) >= 2 and num_heterocyclic_reactants >= 2 and num_complex_reactants >= 2):
                        # Add structural constraint for heterocyclic_reactants_in_coupling_step
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "heterocyclic_reactants_in_coupling_step",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        # Add structural constraint for complex_reactants_in_coupling_step
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants_in_coupling_step (>=8 atoms)",
                                "operator": ">=",
                                "value": 2
                            }
                        })

                        # Verify the product also contains heterocycles
                        product_heterocycles = 0
                        for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, product):
                                product_heterocycles += 1
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                        if product_heterocycles >= 1:
                            # Add structural constraint for heterocyclic_product_in_coupling_step
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "heterocyclic_product_in_coupling_step",
                                    "operator": ">=",
                                    "value": 1
                                }
                            })
                            found_convergent_coupling = True
            except Exception as e:
                # Silently ignore errors for robustness in large-scale processing
                pass

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Remove duplicate entries from atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    # Remove duplicate entries from structural_constraints list (based on content)
    unique_constraints = []
    seen_constraints = set()
    for constraint in findings_json["structural_constraints"]:
        constraint_str = str(constraint) # Convert dict to string for set comparison
        if constraint_str not in seen_constraints:
            unique_constraints.append(constraint)
            seen_constraints.add(constraint_str)
    findings_json["structural_constraints"] = unique_constraints

    return found_convergent_coupling, findings_json
