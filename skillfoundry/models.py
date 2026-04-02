"""Core data models for the codex-exec automation framework."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(slots=True)
class FocusLeaf:
    """A ranked taxonomy leaf that the framework should prioritize."""

    taxonomy_key: str
    domain_slug: str
    domain_name: str
    leaf_name: str
    topic_slug: str
    coverage_status: str
    skill_count: int
    resource_count: int
    verified_skill_count: int
    skill_slugs: list[str] = field(default_factory=list)
    resource_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(slots=True)
class RepositorySummary:
    """Lightweight summary of the current repository state."""

    resource_count: int
    skill_count: int
    taxonomy_domain_count: int
    covered_leaf_count: int
    frontier_leaf_count: int
    todo_leaf_count: int
    covered_domain_count: int
    empty_domain_count: int
    focus_leaves: list[FocusLeaf] = field(default_factory=list)

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["focus_leaves"] = [leaf.to_dict() for leaf in self.focus_leaves]
        return payload


@dataclass(slots=True)
class CodexInvocationResult:
    """Captured result from a single `codex exec` invocation."""

    command: list[str]
    returncode: int
    duration_seconds: float
    stdout: str
    stderr: str
    raw_message: str
    parsed_message: dict
    prompt_path: Path
    schema_path: Path
    stdout_path: Path
    stderr_path: Path
    last_message_path: Path

    def to_dict(self) -> dict:
        return {
            "command": self.command,
            "returncode": self.returncode,
            "duration_seconds": self.duration_seconds,
            "stdout_tail": self.stdout.strip().splitlines()[-20:],
            "stderr_tail": self.stderr.strip().splitlines()[-20:],
            "stdout_path": str(self.stdout_path),
            "stderr_path": str(self.stderr_path),
            "prompt_path": str(self.prompt_path),
            "schema_path": str(self.schema_path),
            "last_message_path": str(self.last_message_path),
            "parsed_message": self.parsed_message,
        }
