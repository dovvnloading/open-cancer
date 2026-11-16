"""Offline CIViC protocol, provenance, permissions, and transport regression tests."""
import copy
from email.message import Message
from email.utils import formatdate
import hashlib
import io
import json
import time
import unittest
from unittest.mock import MagicMock, Mock, patch
from urllib.error import HTTPError, URLError

from scripts import civic_evidence as civic


def api_page(nodes, *, more=False, offset=0, total=None):
    if total is None:
        total = offset + len(nodes) + int(more)
    return json.dumps({"data": {"evidenceItems": {"totalCount": total, "nodes": nodes, "pageInfo": {
        "hasNextPage": more, "endCursor": civic._cursor(offset + len(nodes)) if nodes else None,
    }}}}).encode()


class CivicEvidenceTest(unittest.TestCase):
    def setUp(self):
        self.catalog = json.loads(civic.CATALOG.read_text())
        self.nodes = json.loads(civic.FIXTURE.read_text())["data"]["evidenceItems"]["nodes"]

    def snapshot(self, *pages, limit=4, size=1):
        transport = Mock(side_effect=pages)
        with patch.object(civic, "_fetch_page", transport), patch.object(civic.time, "sleep"):
            bundle = civic.retrieve(self.catalog, page_size=size, max_pages=limit)
        return bundle, transport

    def test_replay_fixture_and_scope(self):
        bundle = civic.retrieve_fixture(self.catalog)
        self.assertIs(bundle, civic.validate_bundle(bundle))
        self.assertEqual(bundle, civic.retrieve_fixture(self.catalog))
        self.assertIs(bundle["synthetic"], True)
        self.assertEqual(len(bundle["records"]), 3)
        job = {"task_id": civic.TASK, "source_ids": ["civic"], "dimensions": {"alteration": "EGFR L858R"}}
        records = civic.records_for_job(bundle, job)
        self.assertEqual([r["evidence_direction"] for r in records], ["SUPPORTS", "DOES_NOT_SUPPORT"])
        for key, value in (("task_id", "unrelated"), ("source_ids", ["civic", "gdc"]), ("source_ids", ["other"]),
                           ("dimensions", {"alteration": "EGFR"}), ("source_ids", ["civic", "other"]),
                           ("dimensions", {"alteration": "EGFR L858R", "gene": "ALK"}),
                           ("dimensions", {"alteration": "EGFR L858R", "disease": "glioma"})):
            self.assertFalse(civic.supports_job({**job, key: value}))
            self.assertEqual(civic.records_for_job(bundle, {**job, key: value}), [])

    def test_pagination_replay_and_provenance(self):
        one = api_page(self.nodes[:1], more=True)
        two = api_page(self.nodes[1:2], offset=1)
        bundle, transport = self.snapshot(one, two)
        self.assertEqual(bundle["status"], "complete")
        self.assertFalse(bundle["synthetic"])
        self.assertEqual([c.args[0] for c in transport.call_args_list],
                         [{"first": 1, "after": None}, {"first": 1, "after": "MQ"}])
        self.assertEqual(bundle["pages"][1]["response_sha256"], hashlib.sha256(two).hexdigest())
        self.assertEqual(bundle["records"][1]["page"], 1)
        record = bundle["records"][1]
        self.assertEqual(record["url"], civic.LINK + str(self.nodes[1]["id"]))
        edited = dict(self.nodes[1], evidenceLevel="A")
        changed, _ = self.snapshot(one, api_page([edited], offset=1))
        self.assertNotEqual(changed["records"][1]["sha256"], record["sha256"])
        self.assertEqual(changed["records"][0]["sha256"], bundle["records"][0]["sha256"])

    def test_limits_partial_and_empty(self):
        partial, _ = self.snapshot(api_page(self.nodes[:1], more=True), limit=1)
        self.assertEqual(partial["status"], "truncated")
        self.assertEqual(len(partial["records"]), 1)
        empty, _ = self.snapshot(api_page([]))
        self.assertEqual((empty["status"], empty["records"]), ("complete", []))
        with patch.object(civic, "_fetch_page") as fetch:
            for kwargs in ({"page_size": 0}, {"page_size": 26}, {"max_pages": 9}, {"max_pages": True}):
                with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                    civic.retrieve(self.catalog, **kwargs)
            fetch.assert_not_called()

    def test_policy_gates_network_and_replay(self):
        bundle = civic.retrieve_fixture(self.catalog)
        for field, val in (("access", "controlled"), ("api_url", civic.ENDPOINT + "?token=x"),
                           ("id", "other"), ("terms_url", "https://example.org"), ("status", "unknown"),
                           ("license_id", "other")):
            catalog = copy.deepcopy(self.catalog)
            (catalog if field in catalog else catalog["reuse"])[field] = val
            with self.subTest(field=field), patch.object(civic, "_fetch_page") as fetch:
                with self.assertRaises(ValueError):
                    civic.retrieve(catalog)
                with self.assertRaises(ValueError):
                    civic.validate_bundle(bundle, catalog)
                fetch.assert_not_called()

    def test_scope_and_no_uncontrolled_fields(self):
        marker = "EXAMPLE_PRIVATE_VALUE"
        row = self.nodes[0]
        extras = {**row, "description": marker, "internalNotes": {"subject": marker}}
        bundle, _ = self.snapshot(api_page([extras]))
        self.assertNotIn(marker, json.dumps(bundle))
        self.assertNotIn("description", bundle["records"][0])
        d, mp = row["disease"], row["molecularProfile"]
        for edit in ({"molecularProfile": {**mp, "name": "EGFR L858R AND EGFR T790M"}},
                     {"molecularProfile": {**mp, "name": marker}},
                     {"molecularProfile": {**mp, "deprecated": True}},
                     {"disease": {**d, "doid": "123"}}, {"disease": {**d, "name": "NSCLC"}},
                     {"disease": None}, {"disease": {**d, "deprecated": True}},
                     {"status": "SUBMITTED"}, {"status": "REJECTED"}, {"flagged": True}):
            with self.subTest(edit=edit):
                bundle, _ = self.snapshot(api_page([{**row, **edit}]))
                self.assertEqual(bundle["records"], [])
                self.assertNotIn(marker, json.dumps(bundle))

    def test_bad_second_page_discards_partial_evidence_and_private_errors(self):
        first = api_page(self.nodes[:1], more=True)
        marker = "EXAMPLE_PRIVATE_VALUE"
        for bad in (b'{"data":null}', b'{"data": {}, "data": {}}', b'[]', b'\xff', b'{}' * civic.MAX_PAGE_BYTES,
                    json.dumps({"errors": [{"message": marker}], "data": {}}).encode(),
                    api_page([{**self.nodes[1], "evidenceType": marker}], offset=1),
                    api_page([{**self.nodes[1], "id": True}], offset=1),
                    api_page(self.nodes[1:]), api_page(self.nodes[1:2], offset=3),
                    api_page([], more=True), b"[" * 5000, civic._Unavailable("rate-limit")):
            with self.subTest(size=type(bad)):
                bundle, _ = self.snapshot(first, bad)
                self.assertEqual((bundle["status"], bundle["records"]), ("unavailable", []))
                self.assertEqual(len(bundle["pages"]), 1)
                self.assertEqual(bundle["failed_variables"], {"first": 1, "after": "MQ"})
                self.assertNotIn(marker, json.dumps(bundle))

    def test_repeated_or_reordered_ids_fail_closed(self):
        for second in (self.nodes[1], self.nodes[0]):
            bundle, _ = self.snapshot(api_page([self.nodes[1]], more=True), api_page([second], offset=1))
            self.assertEqual((bundle["failure"], bundle["records"]), ("source-changed", []))

    def test_count_drift_or_vanishing_last_page_fails_closed(self):
        first = api_page(self.nodes[:1], more=True, total=3)
        changed = api_page(self.nodes[1:2], offset=1, total=2)
        gone = api_page([], offset=1, total=1)
        for second in (changed, gone):
            bundle, _ = self.snapshot(first, second)
            self.assertEqual((bundle["failure"], bundle["records"]), ("source-changed", []))
            self.assertEqual(bundle["failed_variables"], {"first": 1, "after": "MQ"})

    def test_replay_rejects_unsafe_or_inconsistent_late_edits(self):
        for change in (
                lambda b: b.update(query="other"), lambda b: b.update(additional_data="private"),
                lambda b: b.update(retrieved_at="yesterday"), lambda b: b.update(synthetic=0),
                lambda b: b.update(status="unavailable", failure="http-error"),
                lambda b: b["records"][0].update(url=civic.ENDPOINT),
                lambda b: b["records"][0].update(sha256="0" * 64),
                lambda b: b["records"][0].update(page=True),
                lambda b: b["pages"][0]["variables"].update(after="secret"),
                lambda b: b["pages"][0].update(response_sha256="nondigest"),
                lambda b: b["pages"][0].update(node_count=0),
                lambda b: b["records"][1].update(significance="NEW_UPSTREAM_TERM")):
            b = civic.retrieve_fixture(self.catalog)
            change(b)
            b["sha256"] = civic._digest({k: v for k, v in b.items() if k != "sha256"})
            with self.subTest(change=change), self.assertRaises(ValueError):
                civic.validate_bundle(b, self.catalog)
        b = civic.retrieve_fixture(self.catalog)
        b["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            civic.validate_bundle(b, self.catalog)


class TransportTest(unittest.TestCase):
    def response(self, raw=b"{}", content_type="application/json", url=civic.ENDPOINT):
        resp = MagicMock()
        resp.__enter__.return_value = resp
        resp.status = 200
        resp.geturl.return_value = url
        resp.headers = Message()
        resp.headers["Content-Type"] = content_type
        resp.read1.side_effect = io.BytesIO(raw).read
        return resp

    def http_error(self, status, retry_after=None):
        headers = Message()
        if retry_after:
            headers["Retry-After"] = retry_after
        return HTTPError(civic.ENDPOINT, status, "private text", headers, io.BytesIO())

    def fetch(self, *outcomes):
        opener = Mock()
        opener.open.side_effect = list(outcomes)
        self.opener = opener
        with patch.object(civic, "build_opener", return_value=opener), patch.object(civic.time, "sleep") as sleep:
            self.sleep = sleep
            return civic._fetch_page({"first": 5, "after": None})

    def test_transport_public_post_and_retry(self):
        self.assertEqual(self.fetch(self.http_error(429, "3"), self.response()), b"{}")
        self.sleep.assert_called_once_with(3)
        request = self.opener.open.call_args.args[0]
        self.assertEqual(request.full_url, civic.ENDPOINT)
        self.assertEqual(request.get_method(), "POST")
        self.assertNotIn("Authorization", request.headers)
        self.assertEqual(json.loads(request.data), {"query": civic.QUERY, "variables": {"first": 5, "after": None}})
        self.assertEqual(self.opener.open.call_args.kwargs, {"timeout": civic.TIMEOUT})
        self.assertIsNone(civic._NoRedirect().redirect_request(None, None, 307, "", {}, "https://example.org"))

    def test_errors_bounded_redacted(self):
        for outcomes, reason, attempts in (
                ([self.http_error(429, "9999")], "rate-limit", 1),
                ([self.http_error(429, formatdate(time.time() + 600, usegmt=True))], "rate-limit", 1),
                ([self.http_error(503, "bad date")], "http-error", 1),
                ([self.http_error(301)], "http-error", 1), ([self.http_error(401)], "http-error", 1),
                ([self.http_error(404)], "http-error", 1),
                ([self.http_error(502) for _ in range(3)], "http-error", 3),
                ([TimeoutError(), URLError("private text"), TimeoutError()], "network-error", 3),
                ([self.response(url="https://example.org")], "http-error", 1),
                ([self.response(content_type="text/html")], "malformed-response", 1),
                ([self.response(raw=b"a" * (civic.MAX_PAGE_BYTES + 1))], "malformed-response", 1)):
            with self.subTest(reason=reason, attempts=attempts), self.assertRaises(civic._Unavailable) as cm:
                self.fetch(*outcomes)
            self.assertEqual(str(cm.exception), reason)
            self.assertEqual(self.opener.open.call_count, attempts)
