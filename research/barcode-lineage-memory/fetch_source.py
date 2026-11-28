#!/usr/bin/env python3
"""Retrieve and authenticate public barcode inputs in an untracked source cache.

Original retrieval code; no upstream analysis software is used. This script
writes source workbooks and barcode arrays only outside the checkout or under
its ignored .compute directory. No scientific outcomes are calculated here.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import struct
import time
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST_PATH = HERE / 'data/source_manifest.json'
DEFAULT_CACHE = ROOT / '.compute/barcode-source'
EXPECTED_COLUMNS = [f'DT{r}_P{p}' for r in range(1, 5) for p in range(1, 3)]


def source_manifest():
    return json.loads(MANIFEST_PATH.read_text(encoding='utf-8'))


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def checked_cache_dir(path):
    path = Path(path).expanduser().resolve()
    if path.is_relative_to(ROOT) and not path.is_relative_to(ROOT / '.compute'):
        raise ValueError('Source data inside the checkout must remain below the ignored .compute directory')
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_range(url, first, last, total=None):
    """Require a precise range response; bounded retries for transient failures."""
    for attempt in range(3):
        try:
            req = Request(url, headers={'Range': f'bytes={first}-{last}', 'Accept-Encoding': 'identity'})
            with urlopen(req, timeout=45) as response:
                if response.status != 206:
                    raise ValueError(f'Expected HTTP 206, got {response.status}')
                match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                if not match or (int(match[1]), int(match[2])) != (first, last):
                    raise ValueError('Server returned the wrong byte range')
                if total is not None and int(match[3]) != total:
                    raise ValueError('Archive size has changed')
                data = response.read(last - first + 2)
                if len(data) != last - first + 1:
                    raise ValueError('Incomplete or oversized range response')
                return data
        except HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt == 2:
                raise
        except (URLError, TimeoutError):
            if attempt == 2:
                raise
        time.sleep(attempt + 1)
    raise RuntimeError('unreachable range retrieval state')


def verify_workbook(path, manifest=None):
    manifest = source_manifest() if manifest is None else manifest
    expected = manifest['workbook']
    path = Path(path)
    if path.stat().st_size != expected['bytes'] or sha256_file(path) != expected['sha256']:
        raise ValueError('Workbook size or SHA256 mismatch')
    return expected['sha256']


def fetch_workbook(cache_dir=DEFAULT_CACHE, workers=4):
    cache_dir = checked_cache_dir(cache_dir)
    m = source_manifest()
    target = cache_dir / m['workbook']['filename']
    if target.exists():
        verify_workbook(target, m)
        return target
    archive, member = m['archive'], m['archive']['member']
    offset = member['local_header_offset']
    header = get_range(archive['url'], offset, offset + 29, archive['bytes'])
    fields = struct.unpack('<4s5H3I2H', header)
    signature, version, flags, method, mtime, mdate, crc, compressed, rawsize, namesize, extrasize = fields
    if signature != b'PK\x03\x04' or method != 8 or flags & 9:
        raise ValueError('Unsupported or changed ZIP local-header encoding')
    if compressed != member['compressed_bytes'] or rawsize != m['workbook']['bytes'] or crc != int(m['workbook']['crc32_hex'], 16):
        raise ValueError('ZIP member size or CRC metadata differs from pinned source')
    name = get_range(archive['url'], offset + 30, offset + 29 + namesize, archive['bytes']).decode('utf-8')
    if name != member['path']:
        raise ValueError('Wrong ZIP member at pinned offset')
    start = offset + 30 + namesize + extrasize
    ranges = [(p, min(p + 2**20 - 1, start + compressed - 1)) for p in range(start, start + compressed, 2**20)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(lambda r: get_range(archive['url'], *r, archive['bytes']), ranges))
    decompressor = zlib.decompressobj(-15)
    raw = decompressor.decompress(b''.join(parts), rawsize + 1) + decompressor.flush()
    if not decompressor.eof or decompressor.unused_data or len(raw) != rawsize:
        raise ValueError('Incomplete or oversized ZIP member')
    if zlib.crc32(raw) & 0xffffffff != crc or hashlib.sha256(raw).hexdigest() != m['workbook']['sha256']:
        raise ValueError('Extracted workbook failed CRC32/SHA256 authentication')
    temp = target.with_suffix('.xlsx.partial')
    temp.write_bytes(raw)
    temp.replace(target)
    print(json.dumps({'workbook': str(target), 'bytes': len(raw), 'sha256': m['workbook']['sha256']}), flush=True)
    return target


def content_fingerprints(counts, barcodes):
    import numpy as np
    counts, barcodes = np.asarray(counts), np.asarray(barcodes)
    if counts.ndim != 2 or counts.dtype.kind not in 'iu':
        raise ValueError('Canonical counts require a two-dimensional integer array')
    if counts.size and (counts.min() < 0 or counts.max() > np.iinfo(np.int64).max):
        raise ValueError('Counts must fit nonnegative signed int64')
    if barcodes.ndim != 1 or barcodes.dtype.kind != 'U' or len(barcodes) != len(counts):
        raise ValueError('Expected matching Unicode barcode vector')
    return {'counts_shape': list(counts.shape),
            'counts_sha256': hashlib.sha256(counts.astype('<i8', copy=False).tobytes(order='C')).hexdigest(),
            'barcodes_sha256': hashlib.sha256(('\n'.join(barcodes.tolist()) + '\n').encode('ascii')).hexdigest()}


def verify_content_fingerprints(counts, barcodes, line, manifest=None):
    manifest = source_manifest() if manifest is None else manifest
    actual = content_fingerprints(counts, barcodes)
    expected = manifest['cache_fingerprints'][line]
    if any(actual[k] != expected[k] for k in actual):
        raise ValueError(f'Canonical source-content fingerprint mismatch for {line}')
    return actual


def read_authenticated_cache(path, line):
    import numpy as np
    m = source_manifest()
    with np.load(path, allow_pickle=False) as source:
        if str(source['source_sha256'].item()) != m['workbook']['sha256']:
            raise ValueError('Source workbook fingerprint label mismatch')
        if source['columns'].tolist() != EXPECTED_COLUMNS:
            raise ValueError('Unexpected count column ordering')
        counts, barcodes = source['counts'].copy(), source['barcodes'].copy()
    verify_content_fingerprints(counts, barcodes, line, m)
    return counts, barcodes


def prepare_workbook(workbook, cache_dir=DEFAULT_CACHE, lines=('HCTbc', 'SW6bc')):
    import numpy as np
    import openpyxl
    cache_dir = checked_cache_dir(cache_dir)
    m = source_manifest()
    verify_workbook(workbook, m)
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    try:
        for line in lines:
            expected = m['cache_fingerprints'][line]
            n = expected['counts_shape'][0]
            counts = np.empty((n, 8), dtype=np.int64)
            barcodes = np.empty(n, dtype='U30')
            sheet = book[f'{line}_Barcode_Data']
            rows = sheet.iter_rows(values_only=True)
            header = next(rows)
            if list(header[1:9]) != EXPECTED_COLUMNS or len(header) != 9:
                raise ValueError(f'Unexpected workbook count headers in {line}')
            seen = set()
            index = 0
            for row in rows:
                if not any(value is not None for value in row):
                    continue
                if index >= n or len(row) != 9:
                    raise ValueError(f'Unexpected row count/width for {line}')
                barcode = str(row[0])
                if re.fullmatch('[ACGT]{30}', barcode) is None or barcode in seen:
                    raise ValueError('Invalid or duplicate artificial barcode')
                seen.add(barcode)
                values = row[1:9]
                if any(isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0 or not np.isfinite(v) or int(v) != v or v > np.iinfo(np.int64).max for v in values):
                    raise ValueError('Invalid source count')
                counts[index] = values
                barcodes[index] = barcode
                index += 1
            if index != n:
                raise ValueError(f'Expected {n} source rows for {line}, found {index}')
            fp = verify_content_fingerprints(counts, barcodes, line, m)
            target = cache_dir / f'{line}.npz'
            with (cache_dir / f'{line}.npz.partial').open('wb') as f:
                np.savez_compressed(f, counts=counts, barcodes=barcodes,
                                    columns=np.array(EXPECTED_COLUMNS),
                                    source_sha256=np.array(m['workbook']['sha256']))
            (cache_dir / f'{line}.npz.partial').replace(target)
            print(json.dumps({'line': line, 'cache': str(target), **fp}), flush=True)
    finally:
        book.close()


def first_header(url, expected_bytes):
    """Read only a bounded gzip prefix and retain no sequence or quality data."""
    prefix = get_range(url, 0, min(65536, expected_bytes) - 1, expected_bytes)
    decoded = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(prefix, 4096)
    if b'\n' not in decoded:
        raise ValueError('No complete FASTQ header within bounded decoded prefix')
    header = decoded.split(b'\n', 1)[0].decode('ascii')
    match = re.search(r'(?:^@| )([A-Za-z0-9_-]+):(\d+):([A-Za-z0-9_-]+):(\d+):(\d+):(-?\d+):(-?\d+)', header)
    if not match:
        raise ValueError('First FASTQ header has no recognizable Illumina instrument fields')
    return {'url': url, 'compressed_prefix_bytes': len(prefix),
            'compressed_prefix_sha256': hashlib.sha256(prefix).hexdigest(),
            'first_read_header': header,
            'instrument': match[1], 'instrument_run': int(match[2]),
            'flowcell': match[3], 'lane': int(match[4])}


def audit_sequencing(cache_dir=DEFAULT_CACHE, workers=4):
    """Fresh name-based first-header audit, intentionally not a count crosswalk."""
    cache_dir = checked_cache_dir(cache_dir)
    fields = 'study_accession,sample_accession,experiment_accession,run_accession,sample_title,experiment_title,library_strategy,library_source,library_layout,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5'
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urlencode({
        'accession': 'PRJEB74046', 'result': 'read_run', 'fields': fields,
        'format': 'json', 'download': 'false'})
    with urlopen(url, timeout=45) as response:
        metadata = response.read(2**21)
        if response.read(1):
            raise ValueError('Unexpectedly large source metadata response')
    (cache_dir / 'ena_metadata.json').write_bytes(metadata)
    rows = json.loads(metadata)
    matches, missing = [], []
    for line in ('HCTbc', 'SW6bc'):
        for culture in range(1, 5):
            for passage in range(1, 3):
                title = f'QR_{line}_DT{culture}_P{passage}_ce1'
                selected = [r for r in rows if r['sample_title'] == title and r['library_strategy'] == 'AMPLICON']
                if len(selected) != 1:
                    missing.append({'cell_line': line, 'column': f'DT{culture}_P{passage}',
                                    'sample_title_query': title, 'match_count': len(selected)})
                    continue
                row = selected[0]
                matches.append({'cell_line': line, 'column': f'DT{culture}_P{passage}',
                                'mapping_basis': 'exact unique sample-title match, AMPLICON strategy',
                                **row})

    def inspect(row):
        result = dict(row)
        files, sizes = row['fastq_ftp'].split(';'), row['fastq_bytes'].split(';')
        pairs = [(u, int(n)) for u, n in zip(files, sizes, strict=True) if u.endswith('_1.fastq.gz')]
        if len(pairs) != 1:
            raise ValueError('Expected exactly one ENA R1 file')
        file, size = pairs[0]
        try:
            result['ena_first_record'] = {'available': True, **first_header('https://' + file, size)}
        except (HTTPError, URLError, TimeoutError, ValueError) as error:
            result['ena_first_record'] = {'available': False, 'url': 'https://' + file,
                'failure': f'HTTP {error.code}' if isinstance(error, HTTPError) else type(error).__name__}
        print(json.dumps({'column': row['cell_line'] + '_' + row['column'],
                          'available': result['ena_first_record']['available'],
                          'flowcell': result['ena_first_record'].get('flowcell'),
                          'failure': result['ena_first_record'].get('failure')}), flush=True)
        return result

    with ThreadPoolExecutor(max_workers=workers) as pool:
        audited = list(pool.map(inspect, matches))
    original = next(r for r in audited if r['cell_line'] == 'HCTbc' and r['column'] == 'DT1_P1')
    pairs = [(u, int(n)) for u, n in zip(original['submitted_ftp'].split(';'), original['submitted_bytes'].split(';'), strict=True)
             if '_R1_' in u]
    if len(pairs) != 1:
        raise ValueError('Expected exactly one original-submitted R1 crosscheck file')
    file, size = pairs[0]
    try:
        original_check = {'available': True, **first_header('https://' + file, size)}
    except (HTTPError, URLError, TimeoutError, ValueError) as error:
        original_check = {'available': False, 'url': 'https://' + file,
                          'failure': f'HTTP {error.code}' if isinstance(error, HTTPError) else type(error).__name__}
    report = {
        'study_accession': 'PRJEB74046', 'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
        'metadata_url': url, 'metadata_bytes': len(metadata),
        'metadata_sha256': hashlib.sha256(metadata).hexdigest(), 'study_run_count': len(rows),
        'scope': 'Exact unique sample-title matches and the first FASTQ record only. No raw-to-workbook mapping, full-run provenance or contamination rate is established.',
        'read_sequences_or_qualities_retained': False,
        'matched_records': audited, 'unresolved_sample_names': missing,
        'original_submitted_crosscheck': {'run_accession': original['run_accession'],
                                        'sample_title': original['sample_title'], **original_check},
    }
    path = cache_dir / 'sequencing_provenance.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'audit_path': str(path), 'matched': len(audited), 'unresolved': len(missing)}), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache-dir', type=Path, default=DEFAULT_CACHE)
    p.add_argument('--workbook', type=Path)
    p.add_argument('--line', choices=('both','HCTbc','SW6bc'), default='both')
    p.add_argument('--workers', type=int, choices=range(1, 5), default=4)
    p.add_argument('--fetch-only', action='store_true')
    p.add_argument('--audit-sequencing', action='store_true', help='Only audit bounded first-record sequencing provenance')
    args = p.parse_args()
    if args.audit_sequencing:
        audit_sequencing(args.cache_dir, args.workers)
        return
    path = args.workbook or fetch_workbook(args.cache_dir, args.workers)
    verify_workbook(path)
    if not args.fetch_only:
        prepare_workbook(path, args.cache_dir, ('HCTbc','SW6bc') if args.line == 'both' else (args.line,))


if __name__ == '__main__':
    main()
