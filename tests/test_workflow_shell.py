"""Offline syntax and behavior tests for runnable GitHub Actions shell blocks.

CI previously validated Python/Node source but missed a truncated bash quote in
the scheduled versioned reviewer archive verification job.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"


def shell_blocks(path: Path):
    """Extract YAML literal run blocks, without third-party YAML tooling."""
    lines = path.read_text(encoding="utf-8").splitlines()
    for i, row in enumerate(lines):
        match = re.fullmatch(r"(\s*)run:\s*\|\s*", row)
        if not match:
            continue
        level = len(match.group(1))
        block = []
        for following in lines[i + 1:]:
            if following.strip() and len(following) - len(following.lstrip()) <= level:
                break
            block.append(following[level + 2:] if following.strip() else "")
        yield "\n".join(block) + "\n"


class WorkflowShellTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("bash"):
            raise unittest.SkipTest("bash is required for release-workflow shell validation")

    def test_every_inline_bash_block_parses(self):
        checked = 0
        for path in sorted(WORKFLOWS.glob("*.yml")):
            for number, script in enumerate(shell_blocks(path), start=1):
                result = subprocess.run(["bash", "-n"], input=script,
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0,
                                 f"{path.name}: run block {number} has bad bash syntax:\n"
                                 f"{result.stderr}")
                checked += 1
        self.assertGreaterEqual(checked, 10, "Expected to check multiple workflows")

    def test_release_monitor_checksum_is_one_strict_download(self):
        script = list(shell_blocks(WORKFLOWS / "reviewer-release.yml"))[-1]
        self.assertEqual(script.count('gh release download "$RELEASE_TAG"'), 1)
        self.assertEqual(script.count("sha256sum --check -"), 1)
        self.assertIn("'^[a-f0-9]{64}$'", script)
        self.assertIn('test -s egfr-luad-review-kit.zip', script)
        self.assertIn("reviewkit-egfr-luad-v1.*", script)

    def test_monitor_rejects_missing_or_wrong_digest_without_network(self):
        script = list(shell_blocks(WORKFLOWS / "reviewer-release.yml"))[-1]
        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp)
            payload = temp / "payload"
            payload.write_bytes(b"public non-sensitive test archive")
            mock_gh = temp / "gh"
            mock_gh.write_text(
                "#!/usr/bin/env bash\n"
                "set -eu\n"
                "if [[ \"$1\" != release ]]; then exit 15; fi\n"
                "if [[ \"$2\" == view ]]; then echo https://example.invalid/release; exit 0; fi\n"
                "if [[ \"$2\" == download ]]; then\n"
                "  while [[ \"$1\" != --dir ]]; do shift; done\n"
                "  cp \"$FAKE_RELEASE_PAYLOAD\" \"$2/egfr-luad-review-kit.zip\"\n"
                "  exit 0\n"
                "fi\n"
                "exit 16\n", encoding="utf-8")
            mock_gh.chmod(0o755)
            env = dict(os.environ, PATH=str(temp) + os.pathsep + os.environ.get("PATH", ""),
                       RUNNER_TEMP=tmp, FAKE_RELEASE_PAYLOAD=str(payload),
                       RELEASE_TAG="reviewkit-egfr-luad-v1.0.0")
            valid = hashlib.sha256(payload.read_bytes()).hexdigest()
            def exercise(digest: str, tag="reviewkit-egfr-luad-v1.0.0"):
                return subprocess.run(["bash", "-e", "-o", "pipefail", "-c", script],
                                      env=dict(env, EXPECTED_SHA=digest, RELEASE_TAG=tag),
                                      text=True, capture_output=True)

            self.assertNotEqual(exercise("").returncode, 0)
            self.assertNotEqual(exercise("x" * 64).returncode, 0)
            self.assertNotEqual(exercise("0" * 64).returncode, 0)
            self.assertNotEqual(exercise(valid, "unreviewed-tag").returncode, 0)
            good = exercise(valid)
            self.assertEqual(good.returncode, 0, good.stderr)
            self.assertIn("OK", good.stdout)


if __name__ == "__main__":
    unittest.main()
