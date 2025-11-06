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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a nitro group is reduced to an amine
    that participates in a late-stage lactamization to form a medium-sized ring.
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

    # Track if nitro group is present in early steps
    nitro_in_early_steps = False
    # Track if lactamization occurs in late steps
    lactam_formation_in_late_steps = False

    # Track molecules with nitro groups
    nitro_containing_molecules = []
    # Track molecules with amines from nitro reduction
    amine_from_nitro_molecules = []

    def dfs_traverse(node, depth=0):
        nonlocal nitro_in_early_steps, lactam_formation_in_late_steps
        nonlocal nitro_containing_molecules, amine_from_nitro_molecules, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for nitro group in molecules (early steps have higher depth)
            if depth > 2:  # Early in synthesis (depth > 2)
                if checker.check_fg("Nitro group", mol_smiles):
                    nitro_in_early_steps = True
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Nitro group", "position": "early_stage"}})
                    print(f"Found nitro group in early step at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro reduction in any step
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction to amine at depth {depth}: {rsmi}")
                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    # The product of this reaction is the amine we want to track.
                    amine_from_nitro_molecules.append(product_smiles)

                # Check for lactamization in late steps
                if depth <= 2:  # Late in synthesis
                    # Check for lactam formation reactions
                    lactam_reaction_found = checker.check_reaction(
                        "Intramolecular amination (heterocycle formation)", rsmi
                    )
                    if lactam_reaction_found:
                        if "Intramolecular amination (heterocycle formation)" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Intramolecular amination (heterocycle formation)")
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Intramolecular amination (heterocycle formation)", "position": "late_stage"}})

                    # First check if any reactant has an amine that came from nitro reduction
                    amine_reactant = None
                    for r in reactants_smiles:
                        if (checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)):
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            if r in amine_from_nitro_molecules:
                                amine_reactant = r
                                break

                    # Then check if product has a lactam (amide in a ring)
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )
                    if product_has_amide:
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    # Check for medium-sized rings in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    medium_rings = False
                    if product_mol:
                        ring_info = product_mol.GetRingInfo()
                        # Check for rings of size 7-12 (medium-sized)
                        for ring_size in range(7, 13):
                            if any(
                                ring_info.IsAtomInRingOfSize(i, ring_size)
                                for i in range(product_mol.GetNumAtoms())
                            ):
                                medium_rings = True
                                if f"ring_size_{ring_size}" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(f"ring_size_{ring_size}")
                                break

                        # If product has amide and medium-sized rings, and reactant had amine from nitro
                        if (
                            amine_reactant
                            and product_has_amide
                            and medium_rings
                            and lactam_reaction_found
                        ):
                            lactam_formation_in_late_steps = True
                            findings_json["structural_constraints"].append({"type": "sequence", "details": {"first": "Reduction of nitro groups to amines", "second": "Intramolecular amination (heterocycle formation)"}})
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["Intramolecular amination (heterocycle formation)", "Amide", "medium_ring"]}})
                            print(f"Found lactam formation in late step at depth {depth}")

                        # Also check if this is a one-pot nitro reduction and lactamization
                        elif (
                            any(checker.check_fg("Nitro group", r) for r in reactants_smiles)
                            and product_has_amide
                            and medium_rings
                        ):
                            lactam_formation_in_late_steps = True
                            if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["Intramolecular amination (heterocycle formation)", "Amide", "medium_ring"]}})
                            print(
                                f"Found one-pot nitro reduction and lactamization in late step at depth {depth}"
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Final determination based on conditions
    result = nitro_in_early_steps and lactam_formation_in_late_steps
    return result, findings_json
