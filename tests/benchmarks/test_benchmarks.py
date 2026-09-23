"""Parser fidelity, DAG-ness and determinism for the benchmark corpus."""

from __future__ import annotations

from pathlib import Path

import pytest

from bkrobust.benchmarks import acquire, bif, bnjson, dagitty, validate
from bkrobust.benchmarks.describe import describe_network, parse_corpus, parse_file, to_mpdag
from bkrobust.benchmarks.parsed import normalise_names

CACHE = Path("results/axisa3/networks")

# Networks small enough to describe (CPDAG included) inside the test budget.
SMALL = ("asia", "cancer", "child", "insurance", "alarm", "sachs")


@pytest.fixture(scope="module")
def corpus():
    try:
        return acquire.load_cached(CACHE)
    except (FileNotFoundError, RuntimeError) as exc:
        pytest.skip(f"benchmark cache unavailable: {exc}")


@pytest.fixture(scope="module")
def networks(corpus):
    parsed, failures = parse_corpus(corpus)
    assert failures == [], f"unexpected parse failures: {failures}"
    return {n.name: n for n in parsed}


def test_corpus_is_complete(corpus):
    assert len(corpus.by_suffix(".bif.gz")) == 24
    assert len(corpus.by_suffix(".txt")) == 12
    assert len(corpus.by_suffix(".json")) == 4
    assert corpus.sdist_sha256 == acquire.SDIST_SHA256


@pytest.mark.parametrize("name", ["asia", "cancer"])
def test_small_bif_round_trips_against_independent_count(corpus, networks, name):
    """The parser's node/edge counts equal a separate line-by-line count."""
    record = next(f for f in corpus.files if f.name == f"{name}.bif.gz")
    counts = validate.count_bif(record.path)
    parsed = networks[name]
    assert parsed.n_nodes == counts.n_variable_blocks
    assert parsed.n_edges == counts.n_parent_references
    assert counts.n_probability_blocks == counts.n_variable_blocks


def test_every_bif_matches_the_independent_count(corpus, networks):
    """Not just the small ones: all 24 benchmarks agree with the second path."""
    for record in corpus.by_suffix(".bif.gz"):
        counts = validate.count_bif(record.path)
        parsed = networks[record.name[: -len(".bif.gz")]]
        assert parsed.n_nodes == counts.n_variable_blocks, record.name
        assert parsed.n_edges == counts.n_parent_references, record.name


def test_asia_has_its_published_structure(networks):
    """A named check on the smallest network, so a silent regression is visible."""
    asia = networks["asia"]
    assert asia.n_nodes == 8
    assert asia.n_edges == 8
    assert set(asia.nodes) == {
        "asia",
        "bronc",
        "dysp",
        "either",
        "lung",
        "smoke",
        "tub",
        "xray",
    }
    assert ("either", "xray") in asia.edges
    assert ("tub", "either") in asia.edges


def test_every_benchmark_is_a_dag(networks):
    """Every BIF and JSON network converts to a fully oriented acyclic MPDAG."""
    for name, parsed in sorted(networks.items()):
        if parsed.source_format == "dagitty" and parsed.bidirected:
            continue
        assert to_mpdag(parsed).is_dag(), name


def test_bidirected_dagitty_file_is_not_treated_as_a_dag(networks):
    """M-bias uses ``<->``; it must be reported as non-DAG, not coerced."""
    m_bias = networks["M-bias"]
    assert len(m_bias.bidirected) == 2
    description = describe_network(m_bias)
    assert description.is_dag is False
    assert description.status == "not_a_dag"
    assert description.n_undirected_edges is None


def test_parsing_is_deterministic(corpus):
    """Two independent parses of the whole corpus give identical objects."""
    first, failures_a = parse_corpus(corpus)
    second, failures_b = parse_corpus(corpus)
    assert failures_a == failures_b
    assert [n.name for n in first] == [n.name for n in second]
    for a, b in zip(first, second):  # noqa: B905 - equal by construction; strict= is 3.10+
        assert a == b


def test_description_is_deterministic(networks):
    """Descriptive fields (timings excluded) repeat exactly."""
    for name in SMALL:
        a = describe_network(networks[name])
        b = describe_network(networks[name])
        assert a.component_sizes == b.component_sizes
        assert (a.n_undirected_edges, a.n_components, a.max_component_size) == (
            b.n_undirected_edges,
            b.n_components,
            b.max_component_size,
        )


def test_small_descriptions_are_internally_consistent(networks):
    for name in SMALL:
        row = describe_network(networks[name])
        assert row.status == "ok"
        assert row.n_undirected_edges is not None
        assert row.n_components == len(row.component_sizes)
        assert row.n_nodes_in_components == sum(row.component_sizes)
        assert row.n_components_ge3 == sum(1 for s in row.component_sizes if s >= 3)
        assert row.max_component_size == (max(row.component_sizes) if row.component_sizes else 0)
        assert 0.0 <= row.undirected_fraction <= 1.0


def test_timeout_is_recorded_not_swallowed(networks):
    """A zero-tolerance budget must surface as a timeout row, not a crash."""
    row = describe_network(networks["insurance"], budget_s=1e-6)
    assert row.status == "timeout"
    assert row.n_undirected_edges is None
    assert row.seconds_to_cpdag is not None


def test_dagitty_rejects_an_unrecognised_line():
    with pytest.raises(dagitty.DagittyParseError):
        dagitty.parse_dagitty_text("dag {\nA -> B\n!!! junk !!!\n}\n", "x", "x.txt", "0")


def test_dagitty_handles_reversed_arrows():
    parsed = dagitty.parse_dagitty_text('dag {\nA [pos="0,0"]\nB <- A\n}\n', "x", "x.txt", "0")
    assert parsed.edges == (("A", "B"),)


def test_bif_rejects_a_non_bif_file():
    with pytest.raises(ValueError, match="no 'variable' blocks"):
        bif.parse_bif_text("nothing here", "x", "x.bif", "0")


def test_bnjson_rejects_a_wrong_schema():
    with pytest.raises(bnjson.BnJsonParseError):
        bnjson.parse_bnjson_text('{"nodes": ["a"]}', "x", "x.json", "0")


def test_name_normalisation_is_injective_and_reversible():
    mapping = normalise_names(["a b", "a_b", "a b"])
    assert len(set(mapping.values())) == len(mapping)
    inverse = {v: k for k, v in mapping.items()}
    assert len(inverse) == len(mapping)


def test_parse_file_rejects_an_unknown_extension(tmp_path):
    path = tmp_path / "thing.dot"
    path.write_text("digraph {}", encoding="utf-8")
    with pytest.raises(ValueError, match="no parser"):
        parse_file(path, "0")
