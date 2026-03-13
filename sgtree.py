#!/usr/bin/env python
"""SGTree - Species tree construction from marker gene phylogenies.

Thin wrapper for backwards compatibility. The actual code lives in the sgtree/ package.
"""
import warnings


warnings.simplefilter("ignore", SyntaxWarning)


from sgtree.cli import main

main()
