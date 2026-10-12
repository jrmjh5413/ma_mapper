"""The API reference has to agree with the code.

Every signature in `docs/api_reference/` is parsed and compared with the real
one. The reference was written for 0.8.x and had drifted: two functions had
gained parameters, one documented a parameter that does not exist, one
advertised four counting modes and a default value that returned nothing, and
two named functions that could not be imported at all. These tests keep that
from happening again, which matters because the reference is what the journal
means by documenting every optional parameter.
"""

import importlib
import inspect
import logging
import pathlib
import re

import pytest

REFERENCE = pathlib.Path(__file__).resolve().parent.parent / "docs" / "api_reference"
SIGNATURE = re.compile(r"^`(ma_mapper\.[A-Za-z_.]+)\((.*)\)`\s*$", re.M)
# `self` belongs to the MafIO method this package monkey-patches in.
IMPLICIT = {"kwargs", "args", "self"}


def documented_signatures():
    found = []
    for page in sorted(REFERENCE.glob("*.md")):
        for match in SIGNATURE.finditer(page.read_text()):
            found.append((page.name, match.group(1), match.group(2)))
    return found


def test_the_reference_exists_and_documents_the_package():
    found = documented_signatures()
    assert len(found) >= 25, f"only {len(found)} signatures found in {REFERENCE}"


@pytest.mark.parametrize(
    "page,dotted,documented",
    [pytest.param(*row, id=f"{row[0]}:{row[1]}") for row in documented_signatures()],
)
def test_documented_signature_matches_the_code(page, dotted, documented):
    module_path, _, name = dotted.rpartition(".")
    module = importlib.import_module(module_path)
    function = getattr(module, name, None)
    assert function is not None, f"{page} documents {dotted}, which does not exist"

    actual = inspect.signature(function)
    documented_names = [p.split("=")[0].strip() for p in documented.split(",") if p.strip()]
    real_names = [p for p in actual.parameters if p not in IMPLICIT]

    undocumented = [n for n in real_names if n not in documented_names]
    assert not undocumented, f"{page} does not document {undocumented} of {dotted}"

    invented = [n for n in documented_names if n.strip("*") not in actual.parameters]
    assert not invented, f"{page} documents {invented}, which {dotted} does not accept"


@pytest.mark.parametrize(
    "page,dotted,documented",
    [pytest.param(*row, id=f"{row[0]}:{row[1]}") for row in documented_signatures()],
)
def test_documented_defaults_match_the_code(page, dotted, documented):
    module_path, _, name = dotted.rpartition(".")
    actual = inspect.signature(getattr(importlib.import_module(module_path), name))
    for part in documented.split(","):
        if "=" not in part:
            continue
        key, _, value = part.partition("=")
        key, value = key.strip(), value.strip()
        if key not in actual.parameters:
            continue
        real = actual.parameters[key].default
        acceptable = {repr(real), str(real), f"'{real}'", f'"{real}"'}
        if isinstance(real, float):  # 0.50 against 0.5
            acceptable |= {f"{real:.2f}", f"{real:g}"}
        if value in acceptable:
            continue
        # a symbolic default, such as logging.INFO for 20, is clearer in the
        # documentation than the number it resolves to
        try:
            if eval(value, {"logging": logging}, {}) == real:  # noqa: S307
                continue
        except Exception:
            pass
        assert value in acceptable, (
            f"{page}: {dotted} documents {key}={value}, the code has {real!r}"
        )


def test_the_maf_count_modes_in_the_reference_are_implemented():
    """The reference listed 'a', 't', 'c', 'g' and defaulted to 'common'.

    None of the five had a branch in the counting loop, so each returned an
    empty list for every window.
    """
    from ma_mapper import extract_maf

    page = (REFERENCE / "06ma_mapper.extract_maf.md").read_text()
    quoted = set(re.findall(r"'([a-z_]+)':", page))
    # the reference groups several modes per bullet, so also take comma lists
    for line in page.splitlines():
        if line.startswith("- '"):
            quoted |= set(re.findall(r"'([a-z_]+)'", line))
    unimplemented = sorted(q for q in quoted if q not in extract_maf.COUNT_ARGS)
    assert not unimplemented, f"the reference advertises unimplemented modes: {unimplemented}"


def test_no_default_writes_into_the_installed_package():
    """A default output path must not point at the package directory.

    It is read-only in a normal install, and it is not where anyone looks.
    """
    package = pathlib.Path(__file__).resolve().parent.parent / "ma_mapper"
    offenders = [
        path.name
        for path in package.glob("*.py")
        if "os.path.dirname(os.path.abspath(__file__))" in path.read_text()
    ]
    assert not offenders, f"these modules still default to the package directory: {offenders}"
