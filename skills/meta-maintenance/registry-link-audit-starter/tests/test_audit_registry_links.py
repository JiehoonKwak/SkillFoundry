from __future__ import annotations

import importlib.util
import io
import json
import subprocess
import tempfile
import urllib.error
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "meta-maintenance"
    / "registry-link-audit-starter"
    / "scripts"
    / "audit_registry_links.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("audit_registry_links", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class RegistryLinkAuditTests(unittest.TestCase):
    def test_fetch_status_retries_transient_http_error(self) -> None:
        class FakeResponse:
            status = 200

            def __init__(self, url: str) -> None:
                self._url = url
                self.headers = {"Content-Type": "text/html"}

            def geturl(self) -> str:
                return self._url

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb) -> None:
                return None

        transient = urllib.error.HTTPError(
            url="https://example.org",
            code=503,
            msg="Service Unavailable",
            hdrs={"Content-Type": "text/html"},
            fp=io.BytesIO(b""),
        )
        with mock.patch.object(MODULE.urllib.request, "urlopen", side_effect=[transient, FakeResponse("https://example.org")]) as urlopen:
            with mock.patch.object(MODULE.time, "sleep"):
                payload = MODULE.fetch_status("https://example.org", attempts=2)
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["status_code"], 200)
        self.assertEqual(urlopen.call_count, 2)

    def test_fetch_status_returns_urlerror_after_retries(self) -> None:
        with mock.patch.object(
            MODULE.urllib.request,
            "urlopen",
            side_effect=urllib.error.URLError("temporary dns failure"),
        ) as urlopen:
            with mock.patch.object(MODULE.time, "sleep"):
                payload = MODULE.fetch_status("https://example.org", attempts=3)
        self.assertFalse(payload["ok"])
        self.assertIsNone(payload["status_code"])
        self.assertIn("URLError", payload["error"])
        self.assertEqual(urlopen.call_count, 3)

    def test_selected_resources_rejects_unknown_ids(self) -> None:
        with self.assertRaises(SystemExit):
            MODULE.selected_resources(["not-a-real-resource"])

    def test_runtime_link_audit_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "audit.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--resource-id",
                    "matplotlib-docs",
                    "--resource-id",
                    "lychee-docs",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["checked_count"], 2)
            self.assertEqual(payload["failing_count"], 0)


if __name__ == "__main__":
    unittest.main()
