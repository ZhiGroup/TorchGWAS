from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class GWASResult:
    table: list[dict]
    run_metadata: dict
    qc_summary: dict
