from __future__ import annotations

import warnings


# ETE3 still emits import-time SyntaxWarning noise on Python 3.12 from old
# string literals. Suppress those warnings so normal SGTree startup stays clean.
warnings.filterwarnings(
    "ignore",
    category=SyntaxWarning,
    module=r"ete3(\.|$)",
)


__version__ = "2.0.0"
