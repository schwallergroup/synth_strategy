"""Utilities for working with LLMs."""

import ast
import contextlib
import io
import re
from typing import List, Tuple

from synth_strategy.utils import check


def extract_imports(code_str: str) -> str:
    """
    Parses the code string using ast and returns a string with all
    top–level import (or from…import) statements.
    """
    lines = code_str.splitlines()
    import_lines = []
    try:
        module = ast.parse(code_str)
    except Exception as e:
        # If parsing fails return an empty string.
        return ""
    for node in module.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            # Use node.lineno and node.end_lineno (Python 3.8+) to capture the lines.
            start = node.lineno - 1
            end = getattr(node, "end_lineno", node.lineno)
            import_lines.extend(lines[start:end])
    return "\n".join(import_lines)


def extract_functions(code_str: str) -> list:
    """
    Parses the code string and returns a list of tuples (func_name, func_code)
    for each top–level function defined in the code.
    """
    funcs: list = []
    try:
        module = ast.parse(code_str)
    except Exception:
        return funcs
    lines = code_str.splitlines()
    for node in module.body:
        if isinstance(node, ast.FunctionDef):
            # node.lineno is 1-based.
            start = node.lineno - 1
            # node.end_lineno is available in Python 3.8+; fallback to grabbing just the header if not.
            end = getattr(node, "end_lineno", start + 1)
            func_code = "\n".join(lines[start:end])
            funcs.append((node.name, func_code))
    return funcs


def test_generated_code_and_capture(
    route_data: dict, code_str: str, checker
) -> Tuple[bool, str, str]:
    """
    Executes the code snippet in a restricted local namespace,
    expects a callable that takes a single argument (route_data).
    Provides access to a Check object as 'checker' in the candidate function.
    Returns: (passed, stdout, error_message)
    """
    # Create a namespace for exec, and insert the checker under 'checker'
    env = {"checker": checker}
    captured_stdout = ""
    error_message = ""
    passed = False
    # try:
    # Redirect stdout so that prints in the executed code will be captured
    with contextlib.redirect_stdout(io.StringIO()) as f:
        try:
            exec(code_str, env)
            # Attempt to get a function called "test_route"
            candidate_func = env.get("test_route")
            # If not found, search for a callable defined by the user.
            if not callable(candidate_func):
                for name, obj in env.items():
                    if callable(obj) and not name.startswith("__"):
                        candidate_func = obj
                        break
            if candidate_func and callable(candidate_func):
                # Call the candidate function using the provided route_data
                result = candidate_func(route_data)
                passed = bool(result)
            else:
                error_message = (
                    "No valid test function found in executed code."
                )
        except Exception as e:
            error_message = str(e)
            passed = False

        captured_stdout = f.getvalue()
    # except Exception as e:
    #     error_message = str(e)
    return passed, captured_stdout, error_message


def parse_generated_code(response_str: str) -> List[str]:
    """
    Find all Python code segments between <generated_code> and </generated_code>.
    Returns them as a list of strings.
    """
    matches = re.findall(r"<code>(.*?)</code>", response_str, flags=re.DOTALL)
    return matches


def parse_extractable_strategy(response_str: str) -> str:
    """
    Find the synthetic strategy within <synthetic_strategy>...</synthetic_strategy> tags.
    If not found, return an empty string.
    """
    match = re.search(
        r"<extractable_strategy>(.*?)</extractable_strategy>",
        response_str,
        flags=re.DOTALL,
    )
    if match:
        return match.group(1)
    return ""
