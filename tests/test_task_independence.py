"""Policy checks for self-contained runnable scientific tasks."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ACTIVE_ROOTS = (PROJECT_ROOT / "scripts", PROJECT_ROOT / "tasks")
TEXT_SUFFIXES = {".py", ".md", ".toml", ".yaml", ".yml"}
BANNED_REFERENCES = ("notebook", ".ipynb")


def test_active_tasks_do_not_name_exploratory_source_documents() -> None:
    violations: list[str] = []
    for root in ACTIVE_ROOTS:
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
                continue
            content = path.read_text(encoding="utf-8").lower()
            matched = [
                reference
                for reference in BANNED_REFERENCES
                if reference in content
            ]
            if matched:
                relative = path.relative_to(PROJECT_ROOT)
                violations.append(
                    f"{relative}: {', '.join(matched)}"
                )

    assert not violations, (
        "Runnable tasks must be source-document independent:\n"
        + "\n".join(violations)
    )
