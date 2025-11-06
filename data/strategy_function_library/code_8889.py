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


SCAFFOLD_FRAGMENTS_SMARTS = ['c1cnccc1', 'C1CNCCN1']

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a key scaffold is assembled from a predefined set of molecular fragments.
    The function checks if all fragments from the `SCAFFOLD_FRAGMENTS_SMARTS` list are present as substructures in the product molecule,
    and that all fragments were also present across the set of reactant molecules. This identifies the specific reaction step
    where the scaffold is formed. The current fragments are: c1cnccc1, C1CNCCN1.
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

    patterns = [Chem.MolFromSmarts(s) for s in SCAFFOLD_FRAGMENTS_SMARTS]
    early_scaffold_assembly = False

    def dfs_traverse(node, depth=0):
        nonlocal early_scaffold_assembly, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                prod_mol = Chem.MolFromSmiles(product)

                # Check if product contains all fragments from the list
                product_contains_all_fragments = False
                if prod_mol:
                    product_contains_all_fragments = all(prod_mol.HasSubstructMatch(p) for p in patterns)
                    if product_contains_all_fragments:
                        # Add detected ring systems to findings_json
                        if 'c1cnccc1' in SCAFFOLD_FRAGMENTS_SMARTS and 'pyridine' not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                        if 'C1CNCCN1' in SCAFFOLD_FRAGMENTS_SMARTS and 'piperazine' not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                if product_contains_all_fragments:
                    # Check if reactants contained all fragments separately
                    fragments_in_reactants = [False] * len(patterns)
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            for i, p in enumerate(patterns):
                                if react_mol.HasSubstructMatch(p):
                                    fragments_in_reactants[i] = True

                    # If all fragments were found across all reactants, it implies assembly
                    if all(fragments_in_reactants):
                        early_scaffold_assembly = True
                        # Add structural constraint to findings_json
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "pyridine",
                                    "piperazine"
                                ],
                                "condition": "All target fragments must be present across the reactant set and as substructures within the product of a single assembly reaction."
                            }
                        })
                        # Add named reaction to findings_json
                        if "scaffold_assembly" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("scaffold_assembly")

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return early_scaffold_assembly, findings_json
