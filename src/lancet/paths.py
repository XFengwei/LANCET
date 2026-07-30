"""Project-root discovery and policy-aware path construction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib


class ProjectRootNotFound(FileNotFoundError):
    """Raised when a LANCET project root cannot be discovered."""


class DataOutputError(ValueError):
    """Raised when an output path points inside the input-only data tree."""


def _is_lancet_project(directory: Path) -> bool:
    pyproject = directory / "pyproject.toml"
    if not pyproject.is_file():
        return False

    with pyproject.open("rb") as handle:
        metadata = tomllib.load(handle)
    return metadata.get("project", {}).get("name") == "lancet"


def _join_under(base: Path, *parts: str | Path) -> Path:
    base = base.resolve()
    candidate = base
    for part in parts:
        component = Path(part)
        if component.is_absolute():
            raise ValueError(f"Path component must be relative: {component}")
        candidate /= component

    candidate = candidate.resolve()
    if not candidate.is_relative_to(base):
        raise ValueError(f"Path escapes its configured root: {candidate}")
    return candidate


@dataclass(frozen=True)
class ProjectPaths:
    """Canonical paths belonging to one LANCET checkout."""

    root: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", self.root.resolve())

    @classmethod
    def discover(cls, start: str | Path | None = None) -> ProjectPaths:
        """Find the project root by walking upward from *start*."""
        location = Path.cwd() if start is None else Path(start)
        location = location.expanduser().resolve()
        if location.is_file():
            location = location.parent

        for candidate in (location, *location.parents):
            if _is_lancet_project(candidate):
                return cls(candidate)

        raise ProjectRootNotFound(
            f"Could not find a LANCET pyproject.toml above {location}"
        )

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def delivered(self) -> Path:
        return self.data / "delivered"

    @property
    def reduction(self) -> Path:
        return self.data / "reduction"

    @property
    def results(self) -> Path:
        return self.root / "results"

    def delivered_input(self, *parts: str | Path) -> Path:
        """Resolve a stable scientific input under ``data/delivered``."""
        return _join_under(self.delivered, *parts)

    def reduction_input(
        self, *parts: str | Path, declared: bool = False
    ) -> Path:
        """Resolve a declared, read-only intermediate reduction input."""
        if not declared:
            raise PermissionError(
                "Reduction inputs require declared=True and must be read-only"
            )
        return _join_under(self.reduction, *parts)

    def result(self, *parts: str | Path) -> Path:
        """Resolve an output under the project's ``results`` directory."""
        return _join_under(self.results, *parts)

    def validate_output(self, path: str | Path) -> Path:
        """Resolve an output path and reject locations inside ``data``."""
        candidate = Path(path).expanduser()
        if not candidate.is_absolute():
            candidate = self.root / candidate
        candidate = candidate.resolve()

        if not candidate.is_relative_to(self.root):
            raise ValueError(f"Output must remain inside the project: {candidate}")
        if candidate.is_relative_to(self.data.resolve()):
            raise DataOutputError(f"Outputs must not be written inside data/: {candidate}")
        return candidate
