import contextlib
import re

def test_generated_code_and_capture(route_data: dict, code_str: str) -> (bool, str, str):
    """
    Executes the code snippet in a restricted local namespace,
    expects a function named `test_route(route_data_dict) -> bool`.

    Returns:
      (passed, stdout, error_message)
    """
    local_ns = {}
    captured_stdout = ""
    error_message = ""
    passed = False
    try:
        with contextlib.redirect_stdout(io.StringIO()) as f:
            exec(code_str, {}, local_ns)
            test_func = local_ns.get("test_route")
            if callable(test_func):
                result = test_func(route_data)
                passed = bool(result)
            captured_stdout = f.getvalue()
    except Exception as e:
        error_message = str(e)
    return passed, captured_stdout, error_message

def parse_generated_code(response_str: str) -> List[str]:
    """
    Find all Python code segments between <generated_code> and </generated_code>.
    Returns them as a list of strings.
    """
    return re.findall(r"<code>(.*?)</code>", response_str, flags=re.DOTALL)

def parse_computational_analysis(response_str: str) -> str:
    """
    Find the synthetic strategy within <synthetic_strategy>...</synthetic_strategy> tags.
    If not found, return an empty string.
    """
    match = re.search(r"<computational_analysis>(.*?)</computational_analysis>", response_str, flags=re.DOTALL)
    if match:
        return match.group(1)
    return ""
