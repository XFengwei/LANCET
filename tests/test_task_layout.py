"""Every task follows the task pattern, writes only to its own results folder,
and tasks/README.md lists the tasks in an order their inputs allow."""

from pathlib import Path
import re
import tomllib

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TASKS = PROJECT_ROOT / "tasks"
# Local-only coordination task (gitignored); it reads diagnostics but runs nothing.
COORDINATION = {"paper_ii_referee_response"}
TASK_NAMES = sorted(
    path.parent.name for path in TASKS.glob("*/task.toml")
    if path.parent.name not in COORDINATION
)
RESULT_REFERENCE = re.compile(r"^results/([A-Za-z0-9_]+)(/|$)")


def _config(name: str) -> dict:
    return tomllib.loads((TASKS / name / "task.toml").read_text(encoding="utf-8"))


def _strings(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)
    elif isinstance(value, str):
        yield value


def task_dependencies(name: str) -> set[str]:
    """Other tasks whose results this task's configuration reads."""
    found = set()
    for text in _strings(_config(name)):
        match = RESULT_REFERENCE.match(text)
        if match and match.group(1) != name:
            found.add(match.group(1))
    return found


def task_levels() -> dict[str, int]:
    """Level 0 needs no other task; a task's level is one above its deepest input."""
    levels: dict[str, int] = {}

    def level(name: str, chain: tuple[str, ...]) -> int:
        if name in chain:
            raise AssertionError("dependency cycle: " + " -> ".join(chain + (name,)))
        if name not in levels:
            needs = task_dependencies(name)
            levels[name] = 1 + max((level(n, chain + (name,)) for n in needs), default=-1)
        return levels[name]

    for name in TASK_NAMES:
        level(name, ())
    return levels


def test_tasks_exist() -> None:
    assert len(TASK_NAMES) >= 20


@pytest.mark.parametrize("name", TASK_NAMES)
def test_task_has_config_runner_and_readme(name: str) -> None:
    for filename in ("task.toml", "run.py", "README.md"):
        assert (TASKS / name / filename).is_file(), filename


@pytest.mark.parametrize("name", TASK_NAMES)
def test_runner_exposes_config_stage_and_overwrite(name: str) -> None:
    source = (TASKS / name / "run.py").read_text(encoding="utf-8")
    for flag in ('"--config"', '"--stage"', '"--overwrite"'):
        assert flag in source, flag


@pytest.mark.parametrize("name", TASK_NAMES)
def test_task_name_matches_folder(name: str) -> None:
    assert _config(name)["task"]["name"] == name


@pytest.mark.parametrize("name", TASK_NAMES)
def test_outputs_stay_in_the_task_results_folder(name: str) -> None:
    config = _config(name)
    root = f"results/{name}"
    assert config["outputs"]["root"] == root
    for text in _strings(config["outputs"]):
        if text.startswith(("results/", "data/")):
            assert text == root or text.startswith(root + "/"), text


@pytest.mark.parametrize("name", TASK_NAMES)
def test_results_read_from_other_tasks_belong_to_existing_tasks(name: str) -> None:
    for other in task_dependencies(name):
        assert other in TASK_NAMES, f"{name} reads results/{other}/, which no task writes"


def test_task_dependencies_have_no_cycle() -> None:
    task_levels()


def _readme_rows() -> dict[str, tuple[int, set[str]]]:
    rows = {}
    text = (TASKS / "README.md").read_text(encoding="utf-8")
    for line in text.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 3 or not cells[0].isdigit():
            continue
        task = cells[1].strip("`")
        needs = set(re.findall(r"`([a-z0-9_]+)`", cells[2]))
        assert task not in rows, f"{task} listed twice"
        rows[task] = (int(cells[0]), needs)
    return rows


def test_tasks_readme_lists_every_task_once() -> None:
    assert set(_readme_rows()) == set(TASK_NAMES)


def test_tasks_readme_needs_match_the_task_inputs() -> None:
    rows = _readme_rows()
    for name in TASK_NAMES:
        assert rows[name][1] == task_dependencies(name), name


def test_tasks_readme_levels_follow_the_dependencies() -> None:
    rows = _readme_rows()
    levels = task_levels()
    for name in TASK_NAMES:
        assert rows[name][0] == levels[name], name
