"""Source-integrity checks on invented data; no network or biological inputs."""
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('barcode_source_test', ROOT / 'research/barcode-lineage-memory/fetch_source.py')
SOURCE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SOURCE)


class SourceAuthenticationTests(unittest.TestCase):
    def setUp(self):
        self.counts = np.array([[3,0,2,1,5,0,7,8],[0,4,0,9,0,6,0,2]], dtype=np.int64)
        self.barcodes = np.array(['A'*30,'C'*30])
        expected = {'counts_shape':[2,8],
                    'counts_sha256':hashlib.sha256(struct.pack('<16q', *self.counts.reshape(-1).tolist())).hexdigest(),
                    'barcodes_sha256':hashlib.sha256(('A'*30+'\n'+'C'*30+'\n').encode()).hexdigest()}
        self.manifest={'workbook':{'sha256':'f'*64},'cache_fingerprints':{'invented':expected}}

    def test_canonical_content_is_portable_but_authenticates_values(self):
        SOURCE.verify_content_fingerprints(self.counts, self.barcodes, 'invented', self.manifest)
        portable = np.asfortranarray(self.counts.astype('>i8'))
        SOURCE.verify_content_fingerprints(portable, self.barcodes, 'invented', self.manifest)
        modified=self.counts.copy(); modified[0,0]+=1
        with self.assertRaisesRegex(ValueError,'fingerprint'):
            SOURCE.verify_content_fingerprints(modified,self.barcodes,'invented',self.manifest)
        with self.assertRaisesRegex(ValueError,'fingerprint'):
            SOURCE.verify_content_fingerprints(self.counts,self.barcodes[::-1],'invented',self.manifest)
        with self.assertRaises(ValueError):
            SOURCE.verify_content_fingerprints(self.counts.astype(float),self.barcodes,'invented',self.manifest)

    def test_cache_reader_rejects_tampering_despite_preserved_metadata(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(SOURCE,'source_manifest',return_value=self.manifest):
            path=Path(directory)/'invented.npz'
            def save(counts,barcodes,columns=SOURCE.EXPECTED_COLUMNS,label='f'*64):
                np.savez(path,counts=counts,barcodes=barcodes,columns=np.array(columns),source_sha256=np.array(label))
            save(self.counts,self.barcodes)
            counts,barcodes=SOURCE.read_authenticated_cache(path,'invented')
            np.testing.assert_array_equal(counts,self.counts)
            modified=self.counts.copy(); modified[0,0]+=1
            save(modified,self.barcodes)
            with self.assertRaisesRegex(ValueError,'fingerprint'):
                SOURCE.read_authenticated_cache(path,'invented')
            save(self.counts,self.barcodes[::-1])
            with self.assertRaisesRegex(ValueError,'fingerprint'):
                SOURCE.read_authenticated_cache(path,'invented')
            save(self.counts,self.barcodes,SOURCE.EXPECTED_COLUMNS[::-1])
            with self.assertRaisesRegex(ValueError,'ordering'):
                SOURCE.read_authenticated_cache(path,'invented')
            save(self.counts,self.barcodes,label='e'*64)
            with self.assertRaisesRegex(ValueError,'label'):
                SOURCE.read_authenticated_cache(path,'invented')

    def test_source_cache_cannot_enter_tracked_checkout(self):
        with self.assertRaisesRegex(ValueError,'ignored .compute'):
            SOURCE.checked_cache_dir(ROOT / 'research/barcode-source')
        with tempfile.TemporaryDirectory(dir=ROOT.parent) as directory:
            self.assertEqual(Path(directory).resolve(),SOURCE.checked_cache_dir(directory))


if __name__ == '__main__':
    unittest.main()
