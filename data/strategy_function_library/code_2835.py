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


HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrrole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrazole", "isoxazole", "isothiazole", "triazole",
    "tetrazole", "pyrimidine", "pyrazine", "pyridazine", "piperidine",
    "piperazine", "morpholine", "thiomorpholine", "indole", "benzimidazole",
    "benzoxazole", "benzothiazole", "quinoline", "isoquinoline",
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Negishi coupling", "Stille reaction_aryl", "Heck terminal vinyl",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Sonogashira alkyne_aryl halide", "Ullmann condensation",
    "Ullmann-Goldberg Substitution amine", "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol", "Goldberg coupling",
    "Stille reaction_vinyl", "Sonogashira acetylene_aryl halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (depth <= 2) convergent coupling reaction that joins two complex fragments. 
    The reaction must be one of the specified `COUPLING_REACTIONS_OF_INTEREST`. 
    One reactant must contain a heterocycle from `HETEROCYCLES_OF_INTEREST`, and another 
    reactant must contain at least two rings.
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

    # Track if we found the convergent synthesis pattern
    found_convergent_synthesis = False

    # Track the depth during traversal
    max_depth_to_check = 2  # Check reactions at depth 0, 1, and 2

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_synthesis, findings_json

        # Check if this is a reaction node with required metadata
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Extract reaction information
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Checking reaction at depth {depth}: {rsmi}", flush=True)

            # Only process late-stage reactions (depth 0, 1, or 2)
            if depth <= max_depth_to_check:
                # Record positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "coupling_reaction",
                        "position": "late_stage",
                        "max_depth": 2
                    }
                })

                reactants_part = rsmi.split(">")[0]

                reactants = reactants_part.split(".")

                # Check if we have at least 2 reactants (convergent)
                if len(reactants) >= 2:
                    print(f"Found convergent reaction with {len(reactants)} reactants", flush=True)
                    # Record count constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants",
                            "operator": ">=",
                            "value": 2,
                            "scope": "reaction_step"
                        }
                    })

                    # Check for coupling reaction
                    is_coupling = False
                    detected_coupling_reaction = None
                    for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Identified as {reaction_type}", flush=True)
                            is_coupling = True
                            detected_coupling_reaction = reaction_type
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

                    if is_coupling:
                        # Check for heterocycle in one reactant
                        has_heterocycle = False
                        heterocycle_reactant = None
                        heterocycle_found = None

                        for reactant in reactants:
                            for heterocycle in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(heterocycle, reactant):
                                    print(
                                        f"Found heterocycle: {heterocycle} in {reactant}",
                                        flush=True,
                                    )
                                    has_heterocycle = True
                                    heterocycle_reactant = reactant
                                    heterocycle_found = heterocycle
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                    break
                            if has_heterocycle:
                                break

                        # Check for complex ring system in another reactant
                        has_complex_ring = False
                        complex_ring_reactant = None

                        for reactant in reactants:
                            # Skip the heterocycle reactant
                            if reactant == heterocycle_reactant:
                                continue

                            try:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    # Consider a complex ring system as having at least 2 rings
                                    ring_count = reactant_mol.GetRingInfo().NumRings()
                                    if ring_count >= 2:
                                        print(
                                            f"Found complex ring system with {ring_count} rings in {reactant}",
                                            flush=True,
                                        )
                                        has_complex_ring = True
                                        complex_ring_reactant = reactant
                                        # Note: The strategy JSON doesn't specify names for 'complex ring systems' beyond count.
                                        # We'll just note the presence of a reactant with >=2 rings.
                                        break
                            except Exception as e:
                                print(f"Error analyzing complex rings: {e}", flush=True)

                        if has_heterocycle and has_complex_ring:
                            print(
                                f"FOUND CONVERGENT SYNTHESIS WITH HETEROCYCLE: {heterocycle_found}",
                                flush=True,
                            )
                            print(f"Heterocycle reactant: {heterocycle_reactant}", flush=True)
                            print(f"Complex ring reactant: {complex_ring_reactant}", flush=True)
                            found_convergent_synthesis = True
                            # Record co-occurrence constraint if met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "reaction_step_reactants",
                                    "targets": [
                                        "reactant_with_specified_heterocycle",
                                        "reactant_with_at_least_2_rings"
                                    ]
                                }
                            })

        # Continue traversing with modified depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_convergent_synthesis, findings_json
