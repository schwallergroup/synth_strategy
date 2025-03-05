from prompts.code_rewriting import rewriting_prompt

import re

REWRITING_PROMPT = rewriting_prompt

def extract_final_code(llm_response: str) -> str:
    """
    Find the text within <final_code>...</final_code> tags.
    If not found, return an empty string.
    """
    match = re.search(r"<final_code>(.*?)</final_code>", llm_response, flags=re.DOTALL)
    if match:
        return match.group(1)
    return ""


async def call_llm_for_code_rewriting(prompt_text: str) -> str:
    """
    Asynchronously send the prompt_text to the rewriting LLM.
    Returns the entire response (including <final_code>).
    """
    # In a real system, you may choose a different model or route.
    # For demonstration, we'll reuse 'router.acompletion' with a smaller temperature:
    try:
        response = await router.acompletion(
            model="deepseek-r1",  # or a different model dedicated to rewriting
            temperature=0.1,
            messages=[{"role": "user", "content": prompt_text}],
        )
        return response.choices[0].message.content
    except Exception as e:
        print(f"Error calling rewriting LLM: {e}", file=sys.stderr)
        return ""


def run_code_and_test_return_bool(code_str: str) -> bool:
    """
    Executes the code in a restricted local namespace and
    expects a function named `test_route(route_data_dict) -> bool`.
    We'll provide a minimal dictionary as route_data_dict for testing.
    If function returns True, pass; else fail.
    """
    local_ns = {}
    route_data_dict = {"dummy_key": "dummy_value"}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(code_str, {}, local_ns)
            test_func = local_ns.get("test_route")
            if not callable(test_func):
                return False
            result = test_func(route_data_dict)
        return bool(result)
    except Exception:
        return False


async def self_iterate(
    model,
    code_and_description: str,
    test_case_pass: str,
    stdout: str,
    errors: str,
    synthetic_strategy: str,
    rewriting_prompt: str,
    rewriting_limit: int = 3
) -> dict:
    """
    Attempts to improve the generated code until it passes
    the "test_route(route_data_dict)->True" or rewriting_limit is reached.

    Returns a dict with:
       "passed": bool
       "attempts": int
       "final_code": str
    """
    rewriting_count = 0
    passed = False
    final_code = ""

    current_code = code_and_description  # Start with the code snippet as "current_code"

    while rewriting_count < rewriting_limit and not passed:
        rewriting_count += 1

        # Build the rewriting prompt
        prepared_prompt = rewriting_prompt.format(
            CODE_AND_DESCRIPTION=current_code,
            TEST_CASE=test_case_pass,
            stdout=stdout,
            ERRORS=errors,
            SYNTHETIC_STRATEGY=synthetic_strategy
        )

        # Call the rewriting LLM
        llm_response = await call_llm_for_code_rewriting(prepared_prompt)
        # Extract any <final_code> block
        final_code_candidate = extract_final_code(llm_response)

        if not final_code_candidate.strip():
            # No code found; break out
            break

        # Test if new code passes
        code_passed = run_code_and_test_return_bool(final_code_candidate)
        if code_passed:
            passed = True
            final_code = final_code_candidate
        else:
            # If it fails, we feed the new code into the next iteration
            current_code = final_code_candidate

    return {
        "passed": passed,
        "attempts": rewriting_count,
        "final_code": final_code
    }