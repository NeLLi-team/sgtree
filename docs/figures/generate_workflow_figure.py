from __future__ import annotations

import argparse
import json
import os
import textwrap
from copy import deepcopy
from pathlib import Path
from urllib import error, request

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTDIR = Path("docs/figures")
OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"
DEFAULT_OPENROUTER_MODEL = "openai/gpt-4.1-mini"

DEFAULT_SPEC = {
    "title": "SGTree Workflow and Benchmarking Design",
    "subtitle": "Marker-gene species-tree inference with RF-guided duplicate cleanup and synthetic contamination benchmarks",
    "core_heading": "Core Inference Pipeline",
    "core_boxes": [
        {"title": "Inputs", "lines": ["Query proteomes", "Marker HMM set"]},
        {"title": "Search and Parse", "lines": ["pyhmmer search", "Hit filtering", "Marker-count matrix"]},
        {"title": "Per-marker Processing", "lines": ["Sequence extraction", "Profile alignment", "Duplicate removal"]},
        {"title": "Supermatrix Build", "lines": ["trimAl trimming", "Marker concatenation"]},
        {"title": "Outputs", "lines": ["Species tree", "Marker counts", "Rendered figure"]},
    ],
    "selection_heading": "Contamination-aware Marker Selection",
    "selection_boxes": [
        {"title": "Marker Trees", "lines": ["Per-marker FastTree runs", "Reference-assisted option"]},
        {"title": "RF-guided Selection", "lines": ["Copy retained if it minimizes", "species-tree RF distance"]},
        {"title": "Optional Singleton Filter", "lines": ["neighbor", "delta-RF", "backbone", "ensemble"]},
        {"title": "Final Tree", "lines": ["Iterative cleanup", "tree_final.nwk"]},
    ],
    "bridge_label": "marker_selection=yes",
    "benchmark_heading": "Benchmark Harness",
    "benchmark_lines": [
        "Truth panels are generated from clean proteomes, then duplicate/triplicate/replacement contamination is injected.",
        "Benchmark scenarios score normalized RF to the clean truth tree, duplicate removal, native retention, and replacement handling.",
    ],
}


def add_box(
    ax,
    x,
    y,
    w,
    h,
    title,
    lines,
    *,
    fc,
    ec="#1f2933",
    lw=1.8,
    title_size=11,
    body_size=9.0,
    title_width=18,
    body_width=26,
    body_rel_y=0.40,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    if title:
        ax.text(
            x + w / 2,
            y + h - 0.04,
            textwrap.fill(title, width=title_width),
            ha="center",
            va="top",
            fontsize=title_size,
            fontweight="bold",
            color="#111827",
        )
    ax.text(
        x + w / 2,
        y + h * body_rel_y,
        "\n".join(textwrap.fill(line, width=body_width) for line in lines),
        ha="center",
        va="center",
        fontsize=body_size,
        color="#334155",
        linespacing=1.35,
    )


def add_arrow(ax, start, end, *, color="#334155", lw=1.6):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=lw,
            color=color,
            shrinkA=6,
            shrinkB=6,
        )
    )


def add_callout(ax, x, y, w, h, title, lines, *, fc="#ffffff", ec="#1f2933"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.8,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.025,
        y + h - 0.035,
        title,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        color="#111827",
    )
    bullet_lines = [f"• {textwrap.fill(line, width=64)}" for line in lines]
    ax.text(
        x + 0.03,
        y + h - 0.07,
        "\n".join(bullet_lines),
        ha="left",
        va="top",
        fontsize=8.8,
        color="#334155",
        linespacing=1.45,
    )


def _normalize_box_list(boxes: list[dict], expected: int) -> list[dict]:
    normalized = []
    for idx in range(expected):
        box = boxes[idx] if idx < len(boxes) else {}
        fallback_title = DEFAULT_SPEC["core_boxes"][idx]["title"] if expected == 5 else DEFAULT_SPEC["selection_boxes"][idx]["title"]
        title = str(box.get("title", fallback_title)).strip() or fallback_title
        lines = box.get("lines")
        fallback_lines = DEFAULT_SPEC["core_boxes"][idx]["lines"] if expected == 5 else DEFAULT_SPEC["selection_boxes"][idx]["lines"]
        if not isinstance(lines, list) or not lines:
            lines = fallback_lines
        cleaned_lines = [str(line).strip() for line in lines[:4] if str(line).strip()]
        if not cleaned_lines:
            cleaned_lines = fallback_lines
        normalized.append(
            {
                "title": title,
                "lines": cleaned_lines,
            }
        )
    return normalized


def normalize_spec(raw: dict | None) -> dict:
    spec = deepcopy(DEFAULT_SPEC)
    if not raw:
        return spec
    for key in ["title", "subtitle", "core_heading", "selection_heading", "bridge_label", "benchmark_heading"]:
        value = raw.get(key)
        if isinstance(value, str) and value.strip():
            spec[key] = value.strip()
    for key in ["benchmark_lines"]:
        value = raw.get(key)
        if isinstance(value, list) and value:
            cleaned = [str(item).strip() for item in value[:3] if str(item).strip()]
            if cleaned:
                spec[key] = cleaned
    if isinstance(raw.get("core_boxes"), list):
        spec["core_boxes"] = _normalize_box_list(raw["core_boxes"], 5)
    if isinstance(raw.get("selection_boxes"), list):
        spec["selection_boxes"] = _normalize_box_list(raw["selection_boxes"], 4)
    spec["title"] = "SGTree Workflow"
    spec["subtitle"] = "Contamination-aware species-tree inference from marker proteins"
    spec["benchmark_heading"] = "Benchmark Harness"
    spec["benchmark_lines"] = [
        "Build clean truth panels from selected genomes and marker subsets.",
        "Inject duplicate, triplicate, and singleton-replacement contamination; score normalized RF and contaminant removal.",
    ]
    return spec


def _request_openrouter(api_key: str, extra_payload: dict | None = None) -> dict | None:
    system_prompt = (
        "You create compact JSON figure specifications for publication-style workflow diagrams. "
        "Return only content that fits the provided schema. Keep labels short, concrete, and scientifically accurate."
    )
    user_prompt = (
        "Create a workflow figure specification for SGTree. "
        "The figure should cover: per-genome proteome inputs, marker HMM inputs, pyhmmer search, hit filtering, "
        "marker-count matrix generation, per-marker extraction and alignment, duplicate removal, trimAl trimming, "
        "supermatrix concatenation, species-tree inference, RF-guided duplicate cleanup, optional singleton filtering "
        "with neighbor/delta-RF/backbone/ensemble modes, final tree inference, and a benchmark harness that injects "
        "duplicate, triplicate, and replacement contamination into clean truth panels. "
        "Keep the structure aligned to five top-row boxes, four second-row boxes, and a final benchmark-harness callout."
    )

    payload = {
        "model": os.environ.get("OPENROUTER_MODEL", DEFAULT_OPENROUTER_MODEL),
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        "plugins": [{"id": "response-healing"}],
    }
    if extra_payload:
        payload.update(extra_payload)

    req = request.Request(
        OPENROUTER_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": "https://github.com/fredlps/sgtree",
            "X-Title": "SGTree Workflow Figure Planner",
        },
        method="POST",
    )
    try:
        with request.urlopen(req, timeout=60) as response:
            body = json.loads(response.read().decode("utf-8"))
    except (error.URLError, TimeoutError, json.JSONDecodeError):
        return None

    try:
        content = body["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError):
        return None

    if isinstance(content, str):
        try:
            return json.loads(content)
        except json.JSONDecodeError:
            return None
    return content if isinstance(content, dict) else None


def load_openrouter_spec() -> dict | None:
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        return None

    schema_payload = {
        "response_format": {
            "type": "json_schema",
            "json_schema": {
                "name": "workflow_figure_spec",
                "strict": True,
                "schema": {
                    "type": "object",
                    "properties": {
                        "title": {"type": "string"},
                        "subtitle": {"type": "string"},
                        "core_heading": {"type": "string"},
                        "core_boxes": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "title": {"type": "string"},
                                    "lines": {"type": "array", "items": {"type": "string"}},
                                },
                                "required": ["title", "lines"],
                                "additionalProperties": False,
                            },
                            "minItems": 5,
                            "maxItems": 5,
                        },
                        "selection_heading": {"type": "string"},
                        "selection_boxes": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "title": {"type": "string"},
                                    "lines": {"type": "array", "items": {"type": "string"}},
                                },
                                "required": ["title", "lines"],
                                "additionalProperties": False,
                            },
                            "minItems": 4,
                            "maxItems": 4,
                        },
                        "bridge_label": {"type": "string"},
                        "benchmark_heading": {"type": "string"},
                        "benchmark_lines": {"type": "array", "items": {"type": "string"}, "minItems": 1, "maxItems": 3},
                    },
                    "required": [
                        "title",
                        "subtitle",
                        "core_heading",
                        "core_boxes",
                        "selection_heading",
                        "selection_boxes",
                        "bridge_label",
                        "benchmark_heading",
                        "benchmark_lines",
                    ],
                    "additionalProperties": False,
                },
            },
        }
    }
    raw = _request_openrouter(api_key, schema_payload)
    if raw is not None:
        return raw

    fallback_payload = {
        "messages": [
            {
                "role": "system",
                "content": (
                    "Return only valid JSON with keys title, subtitle, core_heading, core_boxes, "
                    "selection_heading, selection_boxes, bridge_label, benchmark_heading, benchmark_lines."
                ),
            },
            {
                "role": "user",
                "content": (
                    "Create the SGTree workflow figure spec as JSON only. "
                    "Use exactly five core_boxes and four selection_boxes. Keep labels short."
                ),
            },
        ]
    }
    return _request_openrouter(api_key, fallback_payload)


def choose_spec(planner: str) -> tuple[dict, str]:
    if planner == "static":
        return deepcopy(DEFAULT_SPEC), "static"
    if planner == "openrouter":
        return normalize_spec(load_openrouter_spec()), "openrouter"
    if planner == "auto":
        raw = load_openrouter_spec()
        if raw is not None:
            return normalize_spec(raw), "openrouter"
        return deepcopy(DEFAULT_SPEC), "static"
    raise ValueError(f"Unsupported planner: {planner}")


def render_workflow_figure(spec: dict, *, planner_used: str) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "workflow_overview_spec.json").write_text(json.dumps({"planner_used": planner_used, "spec": spec}, indent=2) + "\n")

    fig, ax = plt.subplots(figsize=(16, 9))
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.03)

    ax.text(0.04, 0.95, spec["title"], fontsize=22, fontweight="bold", color="#111827", ha="left")
    ax.text(0.04, 0.915, spec["subtitle"], fontsize=11.2, color="#475569", ha="left")

    ax.text(0.04, 0.84, spec["core_heading"], fontsize=13, fontweight="bold", color="#111827")
    core_layout = [
        (0.04, 0.60, 0.15, 0.18, "#ffffff"),
        (0.22, 0.60, 0.16, 0.18, "#ffffff"),
        (0.41, 0.60, 0.18, 0.18, "#ffffff"),
        (0.62, 0.60, 0.16, 0.18, "#ffffff"),
        (0.81, 0.60, 0.15, 0.18, "#ffffff"),
    ]
    for (x, y, w, h, color), box in zip(core_layout, spec["core_boxes"], strict=True):
        add_box(
            ax,
            x,
            y,
            w,
            h,
            box["title"],
            box["lines"],
            fc=color,
            title_size=10.8,
            body_size=8.8,
            title_width=16,
            body_width=22,
            body_rel_y=0.36,
        )
    for idx in range(len(core_layout) - 1):
        x, y, w, h, _ = core_layout[idx]
        nx, ny, nw, nh, _ = core_layout[idx + 1]
        add_arrow(ax, (x + w, y + h / 2), (nx, ny + nh / 2))

    ax.text(0.04, 0.49, spec["selection_heading"], fontsize=13, fontweight="bold", color="#111827")
    selection_layout = [
        (0.06, 0.22, 0.22, 0.19, "#ffffff"),
        (0.32, 0.22, 0.22, 0.19, "#ffffff"),
        (0.58, 0.22, 0.20, 0.19, "#ffffff"),
        (0.82, 0.22, 0.14, 0.19, "#ffffff"),
    ]
    for (x, y, w, h, color), box in zip(selection_layout, spec["selection_boxes"], strict=True):
        add_box(
            ax,
            x,
            y,
            w,
            h,
            box["title"],
            box["lines"],
            fc=color,
            title_size=10.5,
            body_size=8.6,
            title_width=15,
            body_width=20,
            body_rel_y=0.20,
        )
    for idx in range(len(selection_layout) - 1):
        x, y, w, h, _ = selection_layout[idx]
        nx, ny, nw, nh, _ = selection_layout[idx + 1]
        add_arrow(ax, (x + w, y + h / 2), (nx, ny + nh / 2))

    add_arrow(ax, (0.50, 0.60), (0.50, 0.43))
    bridge = FancyBboxPatch(
        (0.43, 0.435),
        0.14,
        0.045,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        linewidth=1.4,
        edgecolor="#1f2933",
        facecolor="#ffffff",
    )
    ax.add_patch(bridge)
    ax.text(0.50, 0.458, spec["bridge_label"], va="center", ha="center", fontsize=9.0, color="#111827")

    add_callout(ax, 0.04, 0.02, 0.92, 0.15, spec["benchmark_heading"], spec["benchmark_lines"])

    fig.savefig(OUTDIR / "workflow_overview.svg", bbox_inches="tight", pad_inches=0.08, facecolor=fig.get_facecolor())
    fig.savefig(OUTDIR / "workflow_overview.png", dpi=300, bbox_inches="tight", pad_inches=0.08, facecolor=fig.get_facecolor())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render the SGTree workflow figure")
    parser.add_argument(
        "--planner",
        choices=["auto", "openrouter", "static"],
        default="auto",
        help="Choose how the figure text spec is produced before local rendering",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spec, planner_used = choose_spec(args.planner)
    render_workflow_figure(spec, planner_used=planner_used)
    print(f"workflow figure rendered with planner={planner_used}")


if __name__ == "__main__":
    main()
