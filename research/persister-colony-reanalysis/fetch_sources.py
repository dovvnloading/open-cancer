#!/usr/bin/env python3
"""Fetch and check public nonclinical research sources; never execute their code.

Standard library only. Large source files must stay outside this repository.
This helper verifies source identity, not the truth or clinical validity of data.
The upstream data retain CC BY 4.0; the repository's code licence does not
relicense them. See data/README.md for attribution and extraction decisions.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile


ARCHIVE_URL = (
    "https://data.mendeley.com/public-files/datasets/mvfm7hs9kw/files/"
    "9ad3d76b-6ef6-4981-a6ca-1809d70f4205/file_downloaded"
)
ARCHIVE_NAME = "Russo_et_al_v2.zip"
ARCHIVE_SHA256 = "2f594e2a76e328d82d954102bc89936edd0932e9d01a8daf8a67e463033d3d4f"
SUPPLEMENT_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/"
    "PMC9279152/supplementaryFiles"
)
# These are hashes of the original member bytes, not newline-normalised text.
SUPPLEMENT_MEMBERS = {
    "41588_2022_1105_MOESM4_ESM.xlsx": (
        102701,
        "12a201b7ef42678daca80c57e972558a383dba9fedf48366a1aa7177cc775b4f",
    ),
    "41588_2022_1105_MOESM1_ESM.pdf": (
        5194367,
        "971e5d468ca2263789400feb02ca73690a197fe88fcf6b143812c0c335a9696c",
    ),
    "41588_2022_1105_Fig3_HTML.jpg": (
        116620,
        "78d45c3a575b0b2bff873547ad2aa83db0f4cb51606cddeec3282c3cde2562a7",
    ),
}
CHUNK = 1024 * 1024
MAX_ARCHIVE_BYTES = 4 * 1024**3
MAX_SUPPLEMENT_BYTES = 128 * 1024**2
MAX_ZIP_ENTRIES = 100000
MAX_SELECTED_BYTES = 16 * 1024**2
RELEVANT_SUFFIXES = {".txt", ".csv", ".tsv", ".xlsx", ".xls", ".ipynb"}


class SourceError(RuntimeError):
    """A source could not be acquired or its pinned identity did not match."""


def sha256_file(path: Path) -> str:
    """Hash exact bytes in bounded memory, including multi-gigabyte archives."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file(path: Path, expected: str) -> str:
    actual = sha256_file(path)
    if actual != expected:
        raise SourceError(
            f"SHA-256 mismatch for {path.name}: expected {expected}; got {actual}"
        )
    return actual


class HTTPSRedirects(urllib.request.HTTPRedirectHandler):
    """Allow ordinary public HTTPS redirects without weakening TLS checks."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        parsed = urllib.parse.urlsplit(newurl)
        if parsed.scheme != "https" or parsed.username or parsed.password:
            raise SourceError("Refused a non-HTTPS or credential-bearing redirect")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(url: str, destination: Path, max_bytes: int, expected: str | None) -> dict:
    """Download atomically; failures never leave a file labelled complete."""
    if destination.exists():
        if expected is None:
            raise SourceError(
                f"{destination} already exists; use --supplement for a local bundle"
            )
        verify_file(destination, expected)
        return {"path": destination.as_posix(), "reused": True}
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.username or parsed.password:
        raise SourceError("Only public HTTPS source URLs are allowed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        opener = urllib.request.build_opener(HTTPSRedirects())
        with opener.open(url, timeout=60) as response:
            if response.status != 200:
                raise SourceError(f"Source returned HTTP {response.status}")
            length = response.headers.get("Content-Length")
            if length and int(length) > max_bytes:
                raise SourceError("Source exceeds the configured download size limit")
            fd, name = tempfile.mkstemp(prefix="source-", suffix=".partial", dir=destination.parent)
            temporary = Path(name)
            digest = hashlib.sha256()
            total = 0
            with os.fdopen(fd, "wb") as output:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    total += len(block)
                    if total > max_bytes:
                        raise SourceError("Source exceeds the configured download size limit")
                    digest.update(block)
                    output.write(block)
            if length and total != int(length):
                raise SourceError("Incomplete download: Content-Length does not match")
            actual = digest.hexdigest()
            if expected is not None and actual != expected:
                raise SourceError(f"Downloaded SHA-256 mismatch: expected {expected}; got {actual}")
            temporary.replace(destination)
            return {
                "path": destination.as_posix(),
                "reused": False,
                "bytes": total,
                "sha256_bytes": actual,
            }
    except urllib.error.HTTPError as error:
        raise SourceError(
            f"Public source returned HTTP {error.code}; no verified source was downloaded. "
            "Do not bypass access controls. An official manually obtained archive can "
            "be checked with --archive."
        ) from error
    except (urllib.error.URLError, TimeoutError) as error:
        raise SourceError(f"Public source download failed: {error}") from error
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def checked_entries(archive: zipfile.ZipFile) -> dict[str, zipfile.ZipInfo]:
    """Reject ambiguous or unsafe paths; never trust names as output paths."""
    items = archive.infolist()
    if len(items) > MAX_ZIP_ENTRIES:
        raise SourceError("Too many ZIP entries")
    result = {}
    for item in items:
        name = item.filename
        path = PurePosixPath(name)
        if (
            not name
            or "\\" in name
            or path.is_absolute()
            or ".." in path.parts
            or any(":" in part for part in path.parts)
            or name in result
        ):
            raise SourceError(f"Unsafe or duplicate ZIP member name: {name!r}")
        mode = item.external_attr >> 16
        if stat.S_ISLNK(mode) or item.flag_bits & 1:
            raise SourceError(f"Encrypted or symlink ZIP entry rejected: {name!r}")
        result[name] = item
    return result


def inspect_archive(path: Path) -> dict:
    """Verify the full V2 archive and list candidate data files without extraction."""
    digest = verify_file(path, ARCHIVE_SHA256)
    with zipfile.ZipFile(path) as archive:
        entries = checked_entries(archive)
        candidates = [
            {
                "member": name,
                "bytes": item.file_size,
                "compressed_bytes": item.compress_size,
                "crc32": f"{item.CRC:08x}",
            }
            for name, item in sorted(entries.items())
            if not item.is_dir() and PurePosixPath(name).suffix.lower() in RELEVANT_SUFFIXES
        ]
    return {
        "source_url": ARCHIVE_URL,
        "source_doi": "10.17632/mvfm7hs9kw.2",
        "licence": "CC BY 4.0",
        "licence_url": "https://creativecommons.org/licenses/by/4.0/",
        "path": path.resolve().as_posix(),
        "bytes": path.stat().st_size,
        "sha256_bytes": digest,
        "candidate_members": candidates,
        "note": "Listed names are candidates, not inferred assays; source code was not executed.",
    }


def extract_supplement(path: Path, output_dir: Path) -> dict:
    """Extract only the three pinned public supplemental files, using basenames."""
    selected = []
    with zipfile.ZipFile(path) as archive:
        entries = checked_entries(archive)
        for name, (expected_size, expected_hash) in SUPPLEMENT_MEMBERS.items():
            if name not in entries:
                raise SourceError(f"Missing exact supplemental member: {name}")
            item = entries[name]
            if item.file_size != expected_size or item.file_size > MAX_SELECTED_BYTES:
                raise SourceError(f"Unexpected supplemental member size: {name}")
            # The member is bounded and CRC-checked by ZipFile.read.
            content = archive.read(item)
            actual = hashlib.sha256(content).hexdigest()
            if actual != expected_hash:
                raise SourceError(f"Supplemental member SHA-256 mismatch: {name}")
            target = output_dir / PurePosixPath(name).name
            if target.exists():
                verify_file(target, expected_hash)
            else:
                with target.open("xb") as output:
                    output.write(content)
            selected.append({"member": name, "bytes": len(content), "sha256_bytes": actual})
    return {
        "source_url": SUPPLEMENT_URL,
        "source_doi": "10.1038/s41588-022-01105-z",
        "path": path.resolve().as_posix(),
        "sha256_bytes": sha256_file(path),
        "selected_members": selected,
        "note": "ZIP packaging may vary; the selected original member bytes are pinned.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="Scratch directory outside this repository")
    archive_group = parser.add_mutually_exclusive_group()
    archive_group.add_argument("--download-archive", action="store_true", help="Download and verify the large Mendeley V2 archive")
    archive_group.add_argument("--archive", type=Path, help="Verify an existing officially obtained V2 ZIP")
    supplement_group = parser.add_mutually_exclusive_group()
    supplement_group.add_argument("--paper-supplement", action="store_true", help="Download the public Europe PMC supplemental bundle")
    supplement_group.add_argument("--supplement", type=Path, help="Use an existing Europe PMC supplemental ZIP")
    args = parser.parse_args(argv)
    if not any((args.download_archive, args.archive, args.paper_supplement, args.supplement)):
        parser.error("Choose --download-archive, --archive, --paper-supplement, or --supplement")
    output_dir = args.output_dir.resolve()
    repository = Path(__file__).resolve().parents[2]
    if output_dir == repository or repository in output_dir.parents:
        parser.error("--output-dir must be outside this repository; do not commit large sources")
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "accessed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scope": "nonclinical research only; no treatment advice",
    }
    try:
        archive_path = args.archive
        if args.download_archive:
            archive_path = output_dir / ARCHIVE_NAME
            download(ARCHIVE_URL, archive_path, MAX_ARCHIVE_BYTES, ARCHIVE_SHA256)
        if archive_path is not None:
            report["mendeley"] = inspect_archive(archive_path)
        supplement_path = args.supplement
        if args.paper_supplement:
            supplement_path = output_dir / "PMC9279152_supplementaryFiles.zip"
            download(SUPPLEMENT_URL, supplement_path, MAX_SUPPLEMENT_BYTES, None)
        if supplement_path is not None:
            report["paper_supplement"] = extract_supplement(supplement_path, output_dir)
        manifest = output_dir / "source_verification.json"
        manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Source byte checks passed; record: {manifest.as_posix()}")
        return 0
    except (SourceError, OSError, ValueError, zipfile.BadZipFile) as error:
        print(f"Source verification failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
