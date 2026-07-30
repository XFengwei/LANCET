"""Basic checks for the LANCET Python package."""

import lancet


def test_package_import() -> None:
    """The package is importable and exposes its initial version."""
    assert lancet.__version__ == "0.1.0"
