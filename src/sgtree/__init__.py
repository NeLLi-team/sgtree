from __future__ import annotations

import warnings


# ETE3 still emits import-time SyntaxWarning noise on Python 3.12 from old
# string literals. Suppress those warnings so normal SGTree startup stays clean.
warnings.filterwarnings(
    "ignore",
    category=SyntaxWarning,
    module=r"ete3(\.|$)",
)
# Some dependency SyntaxWarnings are emitted during source compilation before
# reliable module metadata is attached.
warnings.filterwarnings(
    "ignore",
    message=r"invalid escape sequence .*",
    category=SyntaxWarning,
)


from sgtree._version import __version__
