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
    Detects a synthetic strategy involving ring opening to form a diazo intermediate,
    followed by diazo decomposition and ring closure to form a bicyclic system.
    """
    # Track if we've found each component of the strategy
    has_diazo_intermediate = False
    has_ring_opening = False
    has_ring_closure = False
    has_bicyclic_system = False

    # Track the main synthetic pathway
    main_pathway_molecules = {}
    diazo_depth = None

    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal has_diazo_intermediate, has_ring_opening, has_ring_closure, has_bicyclic_system, diazo_depth

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Store molecule in main pathway for analysis
            if is_main_path and not node.get("in_stock", False):
                main_pathway_molecules[depth] = mol_smiles

                # Check for diazo group using the checker function
                if checker.check_fg("Diazo", mol_smiles):
                    has_diazo_intermediate = True
                    diazo_depth = depth
                    print(f"Found diazo intermediate at depth {depth}: {mol_smiles}")

            # Check for bicyclic system in the final product (depth 0)
            if depth == 0:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    ring_info = mol.GetRingInfo()
                    if ring_info.NumRings() >= 2:
                        # Check if rings share atoms (indicating fused rings)
                        ring_atoms = [set(ring) for ring in ring_info.AtomRings()]
                        for i in range(len(ring_atoms)):
                            for j in range(i + 1, len(ring_atoms)):
                                if ring_atoms[i].intersection(ring_atoms[j]):
                                    has_bicyclic_system = True
                                    print(f"Found bicyclic system in final product: {mol_smiles}")
                                    break
                            if has_bicyclic_system:
                                break

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for ring closure reactions involving diazo compounds
            if (
                checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkyne", rsmi)
                or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkene", rsmi)
                or checker.check_reaction(
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne", rsmi
                )
                or checker.check_reaction(
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkene", rsmi
                )
                or checker.check_reaction("Michael-induced ring closure from diazoalkane", rsmi)
            ):
                has_ring_closure = True
                print(f"Found ring closure reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has a ring that's not in the product
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and any(reactant_mols):
                # Check for ring opening reactions that form diazo compounds
                reactant_rings_total = sum(r.GetRingInfo().NumRings() for r in reactant_mols if r)
                product_rings = product_mol.GetRingInfo().NumRings()

                # Check if this is a ring opening reaction forming a diazo group
                if reactant_rings_total > product_rings:
                    # Check if diazo is formed (present in product but not in reactants)
                    if checker.check_fg("Diazo", product) and not any(
                        checker.check_fg("Diazo", r) for r in reactants
                    ):
                        has_ring_opening = True
                        print(f"Found ring opening to form diazo at depth {depth}: {rsmi}")
                    # Also check for hydrazone oxidation to diazoalkane
                    elif checker.check_reaction("Hydrazone oxidation to diazoalkane", rsmi):
                        has_ring_opening = True
                        print(f"Found hydrazone oxidation to diazoalkane at depth {depth}: {rsmi}")

                # Check if this is a ring closure reaction forming a bicyclic system
                if product_rings > reactant_rings_total:
                    # Check if the product has a bicyclic system
                    ring_info = product_mol.GetRingInfo()
                    ring_atoms = [set(ring) for ring in ring_info.AtomRings()]
                    for i in range(len(ring_atoms)):
                        for j in range(i + 1, len(ring_atoms)):
                            if ring_atoms[i].intersection(ring_atoms[j]):
                                has_ring_closure = True
                                print(
                                    f"Found ring closure forming bicyclic system at depth {depth}: {rsmi}"
                                )
                                break
                        if has_ring_closure:
                            break

        # Traverse children
        for i, child in enumerate(node.get("children", [])):
            # First child is considered part of the main pathway in retrosynthesis
            child_is_main_path = is_main_path and (i == 0)
            dfs_traverse(child, depth + 1, child_is_main_path)

    # Start traversal
    dfs_traverse(route)

    # Check if benzyl ester is maintained throughout the main pathway
    # Only check molecules between diazo intermediate and final product
    has_benzyl_ester_throughout = True
    start_depth = 0
    end_depth = diazo_depth if diazo_depth is not None else max(main_pathway_molecules.keys())

    for depth in range(start_depth, end_depth + 1):
        if depth in main_pathway_molecules:
            mol_smiles = main_pathway_molecules[depth]
            # Check for benzyl ester specifically (not just any ester)
            if not (checker.check_fg("Ester", mol_smiles) and "OCc1ccc" in mol_smiles):
                # Allow for the final product to not have a benzyl ester
                if depth > 0:  # Skip the target molecule
                    has_benzyl_ester_throughout = False
                    print(f"Missing benzyl ester at depth {depth} in main pathway: {mol_smiles}")

    # The strategy is present if all components are found
    strategy_present = (
        has_diazo_intermediate
        and has_ring_opening
        and has_ring_closure
        and has_benzyl_ester_throughout
        and has_bicyclic_system
    )

    print(f"Diazo intermediate ring transformation strategy detected: {strategy_present}")
    print(
        f"Components: diazo={has_diazo_intermediate}, ring_opening={has_ring_opening}, "
        f"ring_closure={has_ring_closure}, benzyl_ester={has_benzyl_ester_throughout}, "
        f"bicyclic={has_bicyclic_system}"
    )

    return strategy_present
