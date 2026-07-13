"""Type stub for the compiled `rna_draw._layout_core` pybind11 extension.

The extension itself is built by `src/layout_core/bindings.cpp`
(CMakeLists.txt) and installed alongside this package; it is not part of
the source tree, so mypy needs this stub to type-check callers.
"""

def plot_coords_turtle(structure: str) -> tuple[list[float], list[float]]: ...
def dump_turtle(structure: str) -> tuple[list[float], list[float]]: ...
def dump_tree(
    structure: str, paired: float = 35.0, unpaired: float = 25.0
) -> list[dict[str, object]]: ...
def dump_detections(
    structure: str, paired: float = 35.0, unpaired: float = 25.0, clearance: float = 1.0
) -> list[tuple[int, int, str]]: ...
def plot_coords_puzzler_resolver_off(structure: str) -> tuple[list[float], list[float]]: ...
def plot_coords_puzzler_sibling_only(structure: str) -> tuple[list[float], list[float]]: ...
def plot_coords_puzzler_sibling_ancestor(structure: str) -> tuple[list[float], list[float]]: ...
def plot_coords_puzzler_full(
    structure: str,
    allow_flipping: bool = False,
    max_config_changes: int = 0,
    clearance: float = 1.0,
) -> tuple[list[float], list[float]]: ...
def dump_change_trace(
    structure: str,
    paired: float = 35.0,
    unpaired: float = 25.0,
    clearance: float = 1.0,
    max_config_changes: int = 25000,
    check_ancestor: bool = False,
) -> list[dict[str, object]]: ...
def plot_coords_puzzler_batch(
    structures: list[str],
    allow_flipping: bool = False,
    max_config_changes: int = 0,
    clearance: float = 1.0,
    num_threads: int = 0,
) -> list[tuple[list[float], list[float], bool]]: ...
def check_overlaps_count(
    x: list[float],
    y: list[float],
    pair_map: list[int],
    node_r: float = 10.0,
    backbone_half_width: float = 7.5,
    pair_half_width: float = 7.5,
    tol: float = 1e-6,
) -> int: ...
def check_overlaps_report(
    x: list[float],
    y: list[float],
    pair_map: list[int],
    node_r: float = 10.0,
    backbone_half_width: float = 7.5,
    pair_half_width: float = 7.5,
    tol: float = 1e-6,
) -> list[tuple[str, str, int, str, int, float, float]]: ...
def check_overlaps_batch(
    xs: list[list[float]],
    ys: list[list[float]],
    pair_maps: list[list[int]],
    node_rs: list[float],
    half_width_factor: float = 0.75,
    tol: float = 1e-6,
    num_threads: int = 0,
) -> list[int]: ...
