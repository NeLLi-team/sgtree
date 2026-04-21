import sys
import warnings


warnings.simplefilter("ignore", SyntaxWarning)

for stream in (sys.stdout, sys.stderr):
    reconfigure = getattr(stream, "reconfigure", None)
    if reconfigure is not None:
        reconfigure(line_buffering=True, write_through=True)


try:
    from sgtree.cli import main
except ModuleNotFoundError as exc:
    if exc.name in {"pyhmmer", "pyrodigal", "hdbscan", "xgboost"}:
        raise SystemExit(
            f"Missing runtime dependency '{exc.name}'. "
            "Run SGTree via `pixi run sgtree ...` or activate the Pixi environment first."
        ) from exc
    raise

main()
