from rdkit import Chem
from typing import Tuple, Dict, List
import copy

# This list can be customized to target different heterocyclic scaffolds.
HETEROCYCLES_OF_INTEREST = ["triazole", "pyrimidine", "pyrazole", "imidazole"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function is a template that detects a synthetic strategy involving sequential
    nucleophilic substitutions on a specific list of heterocyclic scaffolds.
    The heterocycles currently targeted are: triazole, pyrimidine, pyrazole, imidazole.
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

    # Track if we've found evidence of the strategy
    substitution_reactions = []
    heterocycle_scaffolds = []

    def dfs_traverse(node, depth=0):
        nonlocal substitution_reactions, heterocycle_scaffolds, findings_json
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocyclic scaffold in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # Check for triazolopyrimidine or similar fused heterocycles
                has_triazole = checker.check_ring("triazole", product_smiles)
                if has_triazole:
                    findings_json["atomic_checks"]["ring_systems"].append("triazole")
                has_pyrimidine = checker.check_ring("pyrimidine", product_smiles)
                if has_pyrimidine:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                has_pyrazole = checker.check_ring("pyrazole", product_smiles)
                if has_pyrazole:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                has_imidazole = checker.check_ring("imidazole", product_smiles)
                if has_imidazole:
                    findings_json["atomic_checks"]["ring_systems"].append("imidazole")

                # Check for individual heterocycles or fused systems
                has_fused_heterocycle = (
                    has_triazole or has_pyrimidine or has_pyrazole or has_imidazole
                )

                if has_fused_heterocycle:
                    heterocycle_scaffolds.append((product_smiles, depth))
                    # print(f"Found heterocyclic scaffold at depth {depth}: {product_smiles}")

                # Check if heterocycle is formed in this reaction
                heterocycle_formed = False
                if (
                    (
                        not any(checker.check_ring("triazole", r) for r in reactants_smiles)
                        and has_triazole
                    )
                    or (
                        not any(checker.check_ring("pyrimidine", r) for r in reactants_smiles)
                        and has_pyrimidine
                    )
                    or (
                        not any(checker.check_ring("pyrazole", r) for r in reactants_smiles)
                        and has_pyrazole
                    )
                    or (
                        not any(checker.check_ring("imidazole", r) for r in reactants_smiles)
                        and has_imidazole
                    )
                ):
                    heterocycle_formed = True
                    # print(f"Found heterocycle formation at depth {depth}: {product_smiles}")

            # Also check reactants for heterocyclic scaffolds
            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if reactant_mol:
                    has_triazole_r = checker.check_ring("triazole", reactant_smiles)
                    if has_triazole_r:
                        findings_json["atomic_checks"]["ring_systems"].append("triazole")
                    has_pyrimidine_r = checker.check_ring("pyrimidine", reactant_smiles)
                    if has_pyrimidine_r:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                    has_pyrazole_r = checker.check_ring("pyrazole", reactant_smiles)
                    if has_pyrazole_r:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                    has_imidazole_r = checker.check_ring("imidazole", reactant_smiles)
                    if has_imidazole_r:
                        findings_json["atomic_checks"]["ring_systems"].append("imidazole")

                    # Check for individual heterocycles or fused systems
                    has_fused_heterocycle_r = (
                        has_triazole_r or has_pyrimidine_r or has_pyrazole_r or has_imidazole_r
                    )

                    if has_fused_heterocycle_r:
                        heterocycle_scaffolds.append((reactant_smiles, depth))
                        # print(
                        #     f"Found heterocyclic scaffold in reactant at depth {depth}: {reactant_smiles}"
                        # )

            # Check for nucleophilic substitution on heterocycle
            is_nucleophilic_sub = False
            if checker.check_reaction("heteroaromatic_nuc_sub", rsmi):
                is_nucleophilic_sub = True
                findings_json["atomic_checks"]["named_reactions"].append("heteroaromatic_nuc_sub")
            if checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi):
                is_nucleophilic_sub = True
                findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_ortho_nitro")
            if checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi):
                is_nucleophilic_sub = True
                findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_para_nitro")
            if checker.check_reaction("N-arylation_heterocycles", rsmi):
                is_nucleophilic_sub = True
                findings_json["atomic_checks"]["named_reactions"].append("N-arylation_heterocycles")

            if is_nucleophilic_sub:
                # Check if any reactant or product has a heterocyclic scaffold
                has_heterocycle_in_reaction = False
                for smiles, _ in heterocycle_scaffolds:
                    if smiles in reactants_smiles or smiles == product_smiles:
                        has_heterocycle_in_reaction = True
                        break

                # Check for heterocycles in this reaction
                has_heterocycle_in_current = any(
                    checker.check_ring(ring, s)
                    for ring in HETEROCYCLES_OF_INTEREST
                    for s in reactants_smiles + [product_smiles]
                )

                # Only count substitutions on heterocycles
                if has_heterocycle_in_reaction or has_heterocycle_in_current:
                    substitution_reactions.append((rsmi, depth))
                    # print(f"Found nucleophilic substitution at depth {depth}: {rsmi}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Analyze the collected data
    # Sort by depth to check if substitutions are sequential
    substitution_reactions.sort(key=lambda x: x[1])
    heterocycle_scaffolds.sort(key=lambda x: x[1])

    # Check if we have at least 2 substitutions
    has_multiple_substitutions = len(substitution_reactions) >= 2

    # Check if we have a heterocyclic scaffold
    has_heterocycle = len(heterocycle_scaffolds) > 0

    # Check if substitutions are sequential (happen at different depths)
    sequential_substitutions = False
    if len(substitution_reactions) >= 2:
        depths = [depth for _, depth in substitution_reactions]
        # Consider sequential if at different depths or if we have multiple substitutions
        sequential_substitutions = len(set(depths)) >= 2 or len(substitution_reactions) >= 2

    # Return True if we found evidence of the strategy
    strategy_detected = has_heterocycle and has_multiple_substitutions and sequential_substitutions

    # Add structural constraint if detected
    if has_multiple_substitutions:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "nucleophilic substitution on a target heterocycle (triazole, pyrimidine, pyrazole, or imidazole)",
                "operator": ">=",
                "value": 2
            }
        })

    # print(f"Heterocycle elaboration strategy detected: {strategy_detected}")
    # print(
    #     f"- Has heterocycle: {has_heterocycle}"
    # )
    # print(
    #     f"- Multiple substitutions: {has_multiple_substitutions} (count: {len(substitution_reactions)})"
    # )
    # print(f"- Sequential substitutions: {sequential_substitutions}")

    # Remove duplicate entries from atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return strategy_detected, findings_json
