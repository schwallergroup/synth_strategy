#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects a late-stage fragment coupling strategy.
    It looks for the joining of two complex fragments in the second half of the synthesis.
    """
    late_coupling_found = False

    # List of common coupling reaction types
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Buchwald-Hartwig",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Heck terminal vinyl",
        "Heck_non-terminal_vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Ullmann condensation",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
    ]

    # Functional groups that indicate complexity
    complexity_fgs = [
        "Aromatic halide",
        "Boronic acid",
        "Boronic ester",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Ketone",
        "Aldehyde",
        "Phenol",
        "Ether",
        "Sulfonamide",
        "Sulfone",
        "Nitro group",
        "Triflate",
        "Tosylate",
        "Mesylate",
        "Guanidine",
        "Urea",
        "Thiourea",
    ]

    # Common rings that indicate complexity
    complexity_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "indole",
        "benzofuran",
        "benzothiophene",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "purine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling_found

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction SMILES: {rsmi}")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth <= 3)
            if len(reactants) > 1 and depth <= 3:
                print(
                    f"Potential late-stage reaction found at depth {depth} with {len(reactants)} reactants"
                )

                # Check if this is a known coupling reaction
                is_coupling = False
                for rxn in coupling_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Detected coupling reaction: {rxn}")
                        is_coupling = True
                        break

                if not is_coupling:
                    # If not a known coupling, check for characteristic patterns
                    if any(checker.check_fg("Aromatic halide", r) for r in reactants) and any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    ):
                        print("Detected potential Suzuki-like coupling pattern")
                        is_coupling = True
                    elif any(checker.check_fg("Aromatic halide", r) for r in reactants) and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    ):
                        print("Detected potential Buchwald-Hartwig-like coupling pattern")
                        is_coupling = True
                    # Check for guanidine-like coupling patterns
                    elif any("C(=N)" in r for r in reactants) and any("NH" in r for r in reactants):
                        print("Detected potential guanidine coupling pattern")
                        is_coupling = True

                # Even if not a typical coupling reaction, check if it's joining complex fragments
                complex_reactants = []
                reactant_mols = []

                for i, r in enumerate(reactants):
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        reactant_mols.append(mol)

                        # Count rings
                        ring_info = mol.GetRingInfo()
                        ring_count = ring_info.NumRings()

                        # Count functional groups
                        fg_count = sum(1 for fg in complexity_fgs if checker.check_fg(fg, r))

                        # Count complex rings
                        complex_ring_count = sum(
                            1 for ring in complexity_rings if checker.check_ring(ring, r)
                        )

                        complexity_score = ring_count + fg_count + complex_ring_count

                        print(
                            f"Reactant {i} complexity: rings={ring_count}, FGs={fg_count}, complex rings={complex_ring_count}, total score={complexity_score}"
                        )

                        # Lower complexity threshold
                        if (ring_count >= 1 and fg_count >= 1) or complexity_score >= 2:
                            print(f"Reactant {i} is complex: {r}")
                            complex_reactants.append((i, r, mol))

                if len(complex_reactants) >= 2:
                    print(f"Found {len(complex_reactants)} complex reactants")

                    # Verify that the product contains significant portions of both reactants
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        product_atom_count = product_mol.GetNumAtoms()
                        reactant_atom_counts = [m.GetNumAtoms() for m in reactant_mols]

                        print(f"Product atom count: {product_atom_count}")
                        print(f"Reactant atom counts: {reactant_atom_counts}")

                        # Check if product is larger than any single reactant
                        if (
                            product_atom_count > max(reactant_atom_counts) * 0.8
                        ):  # Allow some atom loss
                            # Try to find MCS between reactants and product to verify coupling
                            for i, (idx1, r1, mol1) in enumerate(complex_reactants):
                                for idx2, r2, mol2 in complex_reactants[i + 1 :]:
                                    print(f"Checking coupling between reactants {idx1} and {idx2}")

                                    # Verify that both fragments are present in the product
                                    mcs1 = rdFMCS.FindMCS(
                                        [mol1, product_mol], completeRingsOnly=True
                                    )
                                    mcs2 = rdFMCS.FindMCS(
                                        [mol2, product_mol], completeRingsOnly=True
                                    )

                                    if mcs1.numAtoms > 0 and mcs2.numAtoms > 0:
                                        mcs1_ratio = mcs1.numAtoms / mol1.GetNumAtoms()
                                        mcs2_ratio = mcs2.numAtoms / mol2.GetNumAtoms()

                                        print(f"MCS ratios: {mcs1_ratio:.2f}, {mcs2_ratio:.2f}")

                                        # Lower MCS ratio threshold to 0.5
                                        if mcs1_ratio > 0.5 and mcs2_ratio > 0.5:
                                            print(
                                                f"Late-stage fragment coupling confirmed at depth {depth}"
                                            )
                                            late_coupling_found = True
                                            return  # Exit early once found

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)
            if late_coupling_found:
                return  # Exit early once found

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {late_coupling_found}")

    return late_coupling_found
