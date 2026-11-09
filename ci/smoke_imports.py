"""Import every module in python/ so a broken environment or syntax error fails fast."""
import importlib
import pathlib
import sys

root = pathlib.Path(__file__).resolve().parent.parent / "python"
sys.path.insert(0, str(root))

failed = []
for path in sorted(root.glob("*.py")):
    try:
        importlib.import_module(path.stem)
        print(f"ok   {path.stem}")
    except Exception as exc:  # noqa: BLE001 - report every failure, not just the first
        print(f"FAIL {path.stem}: {type(exc).__name__}: {exc}")
        failed.append(path.stem)

sys.exit(1 if failed else 0)
