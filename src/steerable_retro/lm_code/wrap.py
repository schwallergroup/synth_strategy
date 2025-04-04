from pydantic import BaseModel
from typing import Callable, List, Optional

class LMFunction(BaseModel):
    name: str = ""
    description: str = ""
    func: Callable = None

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)