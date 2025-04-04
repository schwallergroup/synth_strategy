import os
import importlib
from .wrap import LMFunction

# IMPORTANT: Each code_n.py file should define a main() function

# Automatically import all code_n.py files from this directory
def import_all_code_files():
    current_dir = os.path.dirname(__file__)
    for filename in os.listdir(current_dir):
        if filename.startswith("code_"):
            module_name = filename[:-3]  # Remove .py extension
            module = importlib.import_module(f".{module_name}", package="steerable_retro.lm_code")
            # Assuming each module has a main() function
            main_function = getattr(module, 'main', None)
            if main_function:
                globals()[module_name] = LMFunction(
                    name=module_name,
                    description=f"Function from {module_name}.\n{main_function.__doc__}",
                    func=main_function,
                )

# Call the function to import all code files
import_all_code_files()