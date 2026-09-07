"""Acquire the benchmark networks and pin exactly which bytes were used.

The networks live inside the pgmpy source distribution as package data
(``pgmpy/utils/example_models/``). We read that sdist directly rather than
installing pgmpy: recent pgmpy pulls a large dependency tree, and we need only
the graph structure, which is a few hundred kilobytes of text.

The canonical upstream host (``bnlearn.com``, which
``pgmpy.utils.get_example_model`` downloads from) is unreachable from this
machine, and guessed GitHub raw URLs 404. The PyPI sdist is the route that
works, and it has the virtue of being content-addressed: PyPI publishes the
sdist's sha256, so the bytes we parse are verifiable against a third party.

Two hashes are recorded, and both matter:

* the sha256 of the whole sdist, which ties our copy to the published artefact;
* the sha256 of every extracted file, which ties each row of the descriptive
  table to a specific network file.

Without the second, a results table cannot distinguish a parsed network from a
network typed in from memory. That distinction is the reason this module exists.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

#: pgmpy release whose sdist ships the example models we parse.
PGMPY_VERSION = "1.0.0"

#: Direct PyPI URL of that sdist, resolved once from the PyPI JSON API and
#: pinned here so a re-download is reproducible rather than API-dependent.
SDIST_URL = (
    "https://files.pythonhosted.org/packages/44/96/"
    "3d168cf7759d634614e3f7f4b27c8718a9a0715e4f5abf60617c8312f117/pgmpy-1.0.0.tar.gz"
)

#: sha256 that PyPI publishes for that sdist. Verified against our local copy.
SDIST_SHA256 = "aef361e0858bbb1de839c54b940b203170609e1822aff37fc6853e715478255a"

#: Path of the example-model directory inside the sdist.
MEMBER_PREFIX = f"pgmpy-{PGMPY_VERSION}/pgmpy/utils/example_models/"

#: Where the extracted networks and the acquisition manifest are cached.
DEFAULT_CACHE = Path("results/axisa3/networks")

_DOWNLOAD_TIMEOUT_S = 120


@dataclass(frozen=True)
class FileRecord:
    """One extracted network file.

    Args:
        name: Base file name, e.g. ``alarm.bif.gz``.
        path: Absolute path of the extracted copy.
        sha256: Hex digest of the extracted bytes.
        size_bytes: Size of the extracted file.
    """

    name: str
    path: Path
    sha256: str
    size_bytes: int


@dataclass(frozen=True)
class Acquisition:
    """The result of acquiring the benchmark corpus.

    Args:
        sdist_path: Path of the sdist that was read.
        sdist_sha256: Hex digest of the sdist bytes.
        sdist_url: URL the sdist came from (recorded even when a cached copy
            was reused, since that copy's hash is checked against it).
        downloaded: Whether this call fetched the sdist over the network.
        cache_dir: Directory the example models were extracted into.
        files: Extracted files, sorted by name.
    """

    sdist_path: Path
    sdist_sha256: str
    sdist_url: str
    downloaded: bool
    cache_dir: Path
    files: tuple[FileRecord, ...]

    def by_suffix(self, suffix: str) -> tuple[FileRecord, ...]:
        """Extracted files whose name ends with ``suffix``, in sorted order."""
        return tuple(f for f in self.files if f.name.endswith(suffix))


def sha256_file(path: Path) -> str:
    """Hex sha256 of a file, read in chunks.

    Args:
        path: File to hash.

    Returns:
        The lowercase hex digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _candidate_sdists(explicit: Path | None, cache_dir: Path) -> list[Path]:
    """Places a previously fetched sdist may already sit, in search order."""
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    candidates.append(cache_dir / f"pgmpy-{PGMPY_VERSION}.tar.gz")
    return [p for p in candidates if p.is_file()]


def locate_or_download_sdist(
    cache_dir: Path = DEFAULT_CACHE,
    sdist_path: Path | None = None,
    allow_download: bool = True,
) -> tuple[Path, str, bool]:
    """Find a local copy of the pgmpy sdist, downloading it if necessary.

    A local copy is accepted only if its sha256 matches :data:`SDIST_SHA256`;
    a mismatching copy is treated as absent rather than trusted, because the
    whole point of the hash is that we parse the published bytes.

    Args:
        cache_dir: Directory the sdist is cached in.
        sdist_path: An explicit path to try before the cache.
        allow_download: Whether a network fetch is permitted when no valid
            local copy is found.

    Returns:
        ``(path, sha256, downloaded)``.

    Raises:
        RuntimeError: If no valid copy exists and downloading is disallowed,
            or if a freshly downloaded file has the wrong digest.
    """
    cache_dir = Path(cache_dir)
    for candidate in _candidate_sdists(sdist_path, cache_dir):
        digest = sha256_file(candidate)
        if digest == SDIST_SHA256:
            return candidate, digest, False

    if not allow_download:
        raise RuntimeError(
            f"no verified pgmpy-{PGMPY_VERSION} sdist found and downloading is disabled"
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / f"pgmpy-{PGMPY_VERSION}.tar.gz"
    tmp = target.with_suffix(".partial")
    with urllib.request.urlopen(SDIST_URL, timeout=_DOWNLOAD_TIMEOUT_S) as response:
        with tmp.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    digest = sha256_file(tmp)
    if digest != SDIST_SHA256:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"downloaded sdist digest {digest} != expected {SDIST_SHA256}")
    tmp.replace(target)
    return target, digest, True


def _safe_member_name(name: str) -> str:
    """Base name of a tar member, refusing anything path-traversal shaped.

    Args:
        name: The member's name inside the archive.

    Returns:
        The base name.

    Raises:
        ValueError: If the member name escapes its directory.
    """
    base = Path(name).name
    if not base or base in {".", ".."} or "/" in base or "\\" in base:
        raise ValueError(f"refusing suspicious tar member name: {name!r}")
    return base


def extract_example_models(sdist: Path, cache_dir: Path) -> tuple[FileRecord, ...]:
    """Extract ``pgmpy/utils/example_models/`` from the sdist into the cache.

    Extraction is flat (base names only) and refuses any member whose name
    would escape the target directory.

    Args:
        sdist: Path of the verified sdist.
        cache_dir: Directory to extract into; ``example_models`` is created
            beneath it.

    Returns:
        A record per extracted file, sorted by name.

    Raises:
        RuntimeError: If the archive contains no example models.
    """
    out_dir = Path(cache_dir) / "example_models"
    out_dir.mkdir(parents=True, exist_ok=True)

    records: list[FileRecord] = []
    with tarfile.open(sdist, "r:gz") as archive:
        members = sorted(
            (m for m in archive.getmembers() if m.isfile() and m.name.startswith(MEMBER_PREFIX)),
            key=lambda m: m.name,
        )
        if not members:
            raise RuntimeError(f"{sdist} contains no members under {MEMBER_PREFIX}")
        for member in members:
            base = _safe_member_name(member.name)
            source = archive.extractfile(member)
            if source is None:  # pragma: no cover - defensive; isfile() checked
                continue
            destination = out_dir / base
            with source, destination.open("wb") as handle:
                shutil.copyfileobj(source, handle)
            records.append(
                FileRecord(
                    name=base,
                    path=destination.resolve(),
                    sha256=sha256_file(destination),
                    size_bytes=destination.stat().st_size,
                )
            )
    return tuple(sorted(records, key=lambda r: r.name))


def acquire(
    cache_dir: Path = DEFAULT_CACHE,
    sdist_path: Path | None = None,
    allow_download: bool = True,
    write_manifest: bool = True,
) -> Acquisition:
    """Locate/verify the sdist, extract the networks, and hash everything.

    Args:
        cache_dir: Where to cache the sdist, the extracted networks and the
            manifest.
        sdist_path: An explicit sdist path to prefer over the cache.
        allow_download: Whether a network fetch is permitted.
        write_manifest: Whether to write ``acquisition_manifest.json``.

    Returns:
        The :class:`Acquisition` record.
    """
    cache_dir = Path(cache_dir)
    sdist, digest, downloaded = locate_or_download_sdist(
        cache_dir=cache_dir, sdist_path=sdist_path, allow_download=allow_download
    )
    files = extract_example_models(sdist, cache_dir)
    result = Acquisition(
        sdist_path=Path(sdist).resolve(),
        sdist_sha256=digest,
        sdist_url=SDIST_URL,
        downloaded=downloaded,
        cache_dir=cache_dir.resolve(),
        files=files,
    )
    if write_manifest:
        write_acquisition_manifest(result, cache_dir / "acquisition_manifest.json")
    return result


def load_cached(cache_dir: Path = DEFAULT_CACHE) -> Acquisition:
    """Rebuild an :class:`Acquisition` from an existing cache, re-hashing every file.

    Lets the parsers and tests run entirely offline against a cache produced by
    an earlier :func:`acquire`, without needing the sdist present. Hashes are
    recomputed from the files on disk rather than trusted from the manifest, and
    a disagreement is an error: the manifest exists to detect exactly that.

    Args:
        cache_dir: Directory containing ``acquisition_manifest.json`` and
            ``example_models/``.

    Returns:
        The reconstructed acquisition record.

    Raises:
        FileNotFoundError: If the manifest or a file it names is missing.
        RuntimeError: If a file's digest differs from the manifest's.
    """
    cache_dir = Path(cache_dir)
    manifest_path = cache_dir / "acquisition_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"no acquisition manifest at {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    records: list[FileRecord] = []
    for entry in sorted(manifest["files"], key=lambda e: str(e["name"])):
        path = cache_dir / "example_models" / str(entry["name"])
        if not path.is_file():
            raise FileNotFoundError(f"manifest names a missing file: {path}")
        digest = sha256_file(path)
        if digest != entry["sha256"]:
            raise RuntimeError(f"{path.name}: digest {digest} != manifest {entry['sha256']}")
        records.append(
            FileRecord(
                name=path.name,
                path=path.resolve(),
                sha256=digest,
                size_bytes=path.stat().st_size,
            )
        )
    return Acquisition(
        sdist_path=Path(str(manifest["sdist_path"])),
        sdist_sha256=str(manifest["sdist_sha256"]),
        sdist_url=str(manifest["sdist_url"]),
        downloaded=False,
        cache_dir=cache_dir.resolve(),
        files=tuple(records),
    )


def write_acquisition_manifest(acquisition: Acquisition, path: Path) -> None:
    """Write the acquisition record as JSON.

    Args:
        acquisition: The record to serialise.
        path: Destination file; parent directories are created.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "pgmpy_version": PGMPY_VERSION,
        "sdist_url": acquisition.sdist_url,
        "sdist_sha256": acquisition.sdist_sha256,
        "sdist_path": str(acquisition.sdist_path),
        "downloaded_this_run": acquisition.downloaded,
        "member_prefix": MEMBER_PREFIX,
        "n_files": len(acquisition.files),
        "files": [
            {"name": f.name, "sha256": f.sha256, "size_bytes": f.size_bytes}
            for f in acquisition.files
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
