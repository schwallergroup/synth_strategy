code_rewriting_prompt="""You are an expert system designed to analyze and optimize chemical process code. Your task is to carefully examine the provided code, ensure it accurately performs its described function, and align it with the given synthetic strategy. Pay close attention to chemical intuition and subtle differences between basic meanings and their chemical implementations.

First, review the following information:

<code_and_description>
{CODE_AND_DESCRIPTION}
</code_and_description>

<test_case>
{TEST_CASE}
</test_case>

<stdout>
{stdout}
</stdout>

<errors>
{ERRORS}
</errors>

Important note: A list of functional group names and reaction type names is appended to the front of the input string. Make sure to account for this in your analysis and modifications.
There is also a list of common RDKIT functions you can use to analyze the code.
This is not a complete set of all RDKIT functions and classes, you may use additional functions and classes as needed.

You have access to the following classes and functions:
- checker.check_fg(fg_name, query_mol_smiles)
- checker.check_reaction(rxn_type_name, query_rxn_smiles)
- checker.get_fg_smart(fg_name)
- checker.get_reaction_smart(rxn_type_name)
- checker.get_fg_atom_indicies(fg_name, query_mol_smiles)
Please perform a rigorous analysis of the code, following these steps:

1. Analyze the function and chemical process:
   Wrap your analysis inside <code_and_chemical_analysis> tags:
   - Summarize the intended purpose of the function based on its description.
   - Examine the code line by line, explaining what each part does.
   - List out each variable and function parameter, explaining their purpose.
   - Identify any discrepancies between the code's actual behavior and its intended purpose.
   - Break down the chemical process step-by-step, considering reactants, products, and intermediate steps.
   - Analyze the chemical intuition behind the process, explaining why certain steps or checks are necessary.
   - Identify potential edge cases or unusual chemical scenarios that the code should handle.

2. Propose modifications (if necessary):
   Wrap your proposal inside <modification_proposal> tags:
   - If the code doesn't match the description or align with the synthetic strategy, suggest specific changes to address these issues.
   - For each proposed change, write out the current code snippet and the proposed modification side by side.
   - Explain the reasoning behind each proposed modification, considering both the technical and chemical aspects.
   - Use checker.check_fg(name, mol_smiles) and checker.check_reaction(name, rxn_smiles) functions instead of making up SMARTs patterns.
   - For reaction types, always use check_reaction directly, unless it's a functional group interconversion.
   - For functional group interconversions, ensure that your code reliably checks that the FG is actually changed in the reaction.
   - Do not modify the function signature; it must always return a single boolean value.
   - Assume all reaction and molecular smiles present in the synthetic route are valid.

3. Error check:
   Wrap your error analysis inside <error_check> tags:
   - If any errors were provided with the code, determine if the function is the cause.
   - If the function is causing errors, identify the problematic sections and propose fixes.

4. Final review:
   Wrap your final review inside <final_review> tags:
   - Summarize all changes made to the code.
   - Confirm that the modified code accurately performs the described task and is grounded in the synthetic strategy.
   - Double-check that all chemical concepts are correctly implemented.
   - Avoid try-except blocks unless necessary for error handling.
   - Ensure the code makes chemical sense and is neither overly simplistic nor overly complicated.
   - Do not try to capture every case; do not test for more than a couple of reaction types unless the strategy requires it.

After completing your analysis, provide the final version of the code in the following format:

<code>
[Insert the final, potentially modified code here. Include print statements for debugging. Include all necessary imports.]
</code>

Remember:
- Think deeply about chemical intuition throughout your analysis and modification process.
- Ensure that your final code rigorously performs the task it is stated to do, taking into account both the basic meaning of the function description and its specific chemical implementation.
- Do not over-engineer the code to pass the tests at the expense of chemical validity.
- Regarding chemical validity, if you are checking for sequences of FGI's/functional groups, make sure to traverse the synthetic route and check for the presence of these groups in the correct order.
- The function must return a single boolean value.

Example output structure:

<code_and_chemical_analysis>
[Your analysis of the function and chemical process goes here]
</code_and_chemical_analysis>

<modification_proposal>
[Your proposed modifications, if any, go here]
</modification_proposal>

<error_check>
[Your error analysis goes here]
</error_check>

<final_review>
[Your final review goes here]
</final_review>

<code>
def example_function(parameter1, parameter2):
    # Your modified code here
    return True  # or False
</code>
"""