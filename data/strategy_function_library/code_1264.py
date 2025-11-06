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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving heterocycle formation (imidazole and pyrazole)
    combined with nitro group chemistry (introduction and reduction).
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

    # Track if we found the key elements of the strategy
    has_nitro_introduction = False
    has_nitro_reduction = False
    has_imidazole_formation = False
    has_pyrazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_introduction, has_nitro_reduction, has_imidazole_formation, has_pyrazole_formation, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check for nitro introduction
                        nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                        reactant_nitro_count = sum(
                            len(r.GetSubstructMatches(nitro_pattern)) for r in reactants if r
                        )
                        product_nitro_count = len(product.GetSubstructMatches(nitro_pattern))

                        if product_nitro_count > reactant_nitro_count:
                            has_nitro_introduction = True
                            findings_json["atomic_checks"]["named_reactions"].append("nitro_group_introduction")
                            findings_json["atomic_checks"]["functional_groups"].append("nitro")
                            # print("Detected nitro group introduction")

                        # Check for nitro reduction
                        amine_pattern = Chem.MolFromSmarts("[NH2]")
                        reactant_amine_count = sum(
                            len(r.GetSubstructMatches(amine_pattern)) for r in reactants if r
                        )
                        product_amine_count = len(product.GetSubstructMatches(amine_pattern))

                        if product_amine_count > reactant_amine_count and product_nitro_count < reactant_nitro_count:
                            has_nitro_reduction = True
                            findings_json["atomic_checks"]["named_reactions"].append("nitro_group_reduction")
                            findings_json["atomic_checks"]["functional_groups"].append("primary_amine")
                            # print("Detected nitro reduction to amine")

                        # Check for imidazole formation
                        imidazole_pattern = Chem.MolFromSmarts("[nH]1cncc1")
                        reactant_imidazole_count = sum(
                            len(r.GetSubstructMatches(imidazole_pattern)) for r in reactants if r
                        )
                        product_imidazole_count = len(
                            product.GetSubstructMatches(imidazole_pattern)
                        )

                        if product_imidazole_count > reactant_imidazole_count:
                            has_imidazole_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append("imidazole_formation")
                            findings_json["atomic_checks"]["ring_systems"].append("imidazole")
                            # print("Detected imidazole formation")

                        # Check for pyrazole formation
                        pyrazole_pattern = Chem.MolFromSmarts("[nH]1ncc[c]1")
                        reactant_pyrazole_count = sum(
                            len(r.GetSubstructMatches(pyrazole_pattern)) for r in reactants if r
                        )
                        product_pyrazole_count = len(product.GetSubstructMatches(pyrazole_pattern))

                        if product_pyrazole_count > reactant_pyrazole_count:
                            has_pyrazole_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append("pyrazole_formation")
                            findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                            # print("Detected pyrazole formation")

                except Exception as e:
                    # print(f"Error processing reaction: {e}")
                    pass # Suppress print for refactored code

        # Process children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # The strategy is present if we have heterocycle formation and nitro chemistry
    strategy_present = (has_imidazole_formation or has_pyrazole_formation) and (
        has_nitro_introduction or has_nitro_reduction
    )

    if strategy_present:
        # print("Detected heterocycle formation with nitro chemistry strategy")
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain at least one heterocycle formation event (imidazole or pyrazole) AND at least one nitro chemistry event (introduction or reduction).",
                "expression": {
                    "operator": "AND",
                    "operands": [
                        {
                            "operator": "OR",
                            "operands": [
                                "imidazole_formation",
                                "pyrazole_formation"
                            ]
                        },
                        {
                            "operator": "OR",
                            "operands": [
                                "nitro_group_introduction",
                                "nitro_group_reduction"
                            ]
                        }
                    ]
                }
            }
        })

    return strategy_present, findings_json