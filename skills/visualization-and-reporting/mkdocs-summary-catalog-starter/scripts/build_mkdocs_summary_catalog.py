#!/usr/bin/env python3
"""Build a tiny MkDocs catalog site and summarize the result."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MKDOCS_BIN = ROOT / "slurm" / "envs" / "reporting" / "bin" / "mkdocs"


def load_catalog(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_workspace(workspace: Path, catalog: dict) -> dict:
    if workspace.exists():
        shutil.rmtree(workspace)
    docs_dir = workspace / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    nav = [{"Home": "index.md"}]
    (docs_dir / "index.md").write_text(f"# {catalog['site_name']}\n\nGenerated starter catalog.\n", encoding="utf-8")
    for page in catalog["pages"]:
        filename = f"{page['slug']}.md"
        nav.append({page["title"]: filename})
        (docs_dir / filename).write_text(f"# {page['title']}\n\n{page['body']}\n", encoding="utf-8")
    mkdocs_yml = workspace / "mkdocs.yml"
    mkdocs_yml.write_text(
        "site_name: " + catalog["site_name"] + "\n" +
        "nav:\n" +
        "".join(
            f"  - {list(item.keys())[0]}: {list(item.values())[0]}\n"
            for item in nav
        ),
        encoding="utf-8",
    )
    return {"nav": nav, "page_count": len(nav)}


def run_build(input_path: Path, workspace: Path) -> dict:
    catalog = load_catalog(input_path)
    workspace_info = write_workspace(workspace, catalog)
    site_dir = workspace / "site"
    completed = subprocess.run(
        [str(MKDOCS_BIN), "build", "--clean", "--site-dir", str(site_dir)],
        cwd=workspace,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    html_files = sorted(path.relative_to(site_dir).as_posix() for path in site_dir.rglob("*.html"))
    return {
        "input_path": str(input_path.resolve()),
        "workspace": str(workspace.resolve()),
        "site_dir": str(site_dir.resolve()),
        "page_count": workspace_info["page_count"],
        "html_file_count": len(html_files),
        "html_files": html_files,
        "stdout_tail": completed.stdout.splitlines()[-3:],
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_build(args.input, args.workspace)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
