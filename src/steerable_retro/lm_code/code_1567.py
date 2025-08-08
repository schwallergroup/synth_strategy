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
    Detects a synthetic strategy involving a late-stage Suzuki coupling
    with convergent fragment assembly.
    """
    suzuki_coupling_found = False
    borylation_found = False
    snar_count = 0
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, borylation_found, snar_count, nitro_reduction_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for Suzuki coupling (depth 0, 1, or 2 = late stage)
                if depth <= 2:
                    if (
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                        or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                        or checker.check_reaction("{Suzuki}", rsmi)
                        or checker.check_reaction("Suzuki", rsmi)
                    ):
                        suzuki_coupling_found = True
                        print(f"Found late-stage Suzuki coupling at depth {depth}")

                # Check for borylation
                if (
                    checker.check_reaction("Preparation of boronic acids", rsmi)
                    or checker.check_reaction("Preparation of boronic ethers", rsmi)
                    or checker.check_reaction(
                        "Preparation of boronic acids from trifluoroborates", rsmi
                    )
                    or checker.check_reaction(
                        "Preparation of boronic acids without boronic ether", rsmi
                    )
                    or checker.check_reaction("Synthesis of boronic acids", rsmi)
                ):
                    borylation_found = True
                    print("Found borylation reaction")

                # Check for SNAr (nucleophilic aromatic substitution)
                if (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                ):
                    snar_count += 1
                    print(f"Found SNAr or Ullmann-type reaction (count: {snar_count})")

                # Additional check for N-arylation which can be similar to SNAr
                if (
                    checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    or checker.check_reaction("{N-arylation_heterocycles}", rsmi)
                    or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                ):
                    snar_count += 1
                    print(f"Found N-arylation reaction (count: {snar_count})")

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_reduction_found = True
                    print("Found nitro reduction")

                # Fallback checks if specific reaction types aren't detected
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Fallback for Suzuki coupling at late stage
                if depth <= 2 and not suzuki_coupling_found:
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    )
                    has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                    if has_boronic and has_aryl_halide:
                        # Check if product has a new C-C bond between aromatic rings
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # More comprehensive biaryl patterns
                            biaryl_pattern1 = Chem.MolFromSmarts("c:c-c:c")
                            biaryl_pattern2 = Chem.MolFromSmarts("c-!:c:c")
                            biaryl_pattern3 = Chem.MolFromSmarts("c:c-c")

                            if (
                                product_mol.HasSubstructMatch(biaryl_pattern1)
                                or product_mol.HasSubstructMatch(biaryl_pattern2)
                                or product_mol.HasSubstructMatch(biaryl_pattern3)
                            ):
                                suzuki_coupling_found = True
                                print(
                                    f"Found late-stage Suzuki coupling (fallback) at depth {depth}"
                                )

                # Fallback for nitro reduction
                if not nitro_reduction_found:
                    has_nitro = any(checker.check_fg("Nitro group", r) for r in reactants)
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Aniline", product)
                        or checker.check_fg("Secondary amine", product)
                    )

                    if has_nitro and has_amine_product:
                        nitro_reduction_found = True
                        print("Found nitro reduction (fallback)")

                # Additional check for Suzuki-like C-C bond formation
                if depth <= 2 and not suzuki_coupling_found:
                    # Check if any reactant has a boronic acid/ester and if the product has a biaryl motif
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and reactant_mols:
                        # Check for biaryl formation
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                        if product_mol.HasSubstructMatch(biaryl_pattern):
                            # Check if this is likely a Suzuki coupling
                            if any(
                                checker.check_fg("Boronic acid", r)
                                or checker.check_fg("Boronic ester", r)
                                for r in reactants
                            ):
                                suzuki_coupling_found = True
                                print(f"Found late-stage Suzuki-like coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # Consider either a Suzuki coupling alone or a combination of borylation with other components
    strategy_present = suzuki_coupling_found or (
        borylation_found and (snar_count >= 1 or nitro_reduction_found)
    )

    print(
        f"Strategy components found: Suzuki={suzuki_coupling_found}, Borylation={borylation_found}, SNAr count={snar_count}, Nitro reduction={nitro_reduction_found}"
    )
    print(f"Strategy present: {strategy_present}")

    return strategy_present
