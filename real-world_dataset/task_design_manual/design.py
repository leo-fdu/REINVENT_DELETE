"""Validate manual LibINVENT/LinkINVENT cuts and write immutable snapshots."""

from __future__ import annotations

from datetime import datetime
from functools import lru_cache
from html import escape
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json
import sys

import torch
from rdkit import Chem, rdBase
from rdkit.Chem.Draw import rdMolDraw2D

HERE = Path(__file__).resolve().parent
ATLAS = HERE.parent / "task_design_astra"
sys.path.insert(0, str(ATLAS))
import render_designs as atlas  # noqa: E402

MODES = ("libinvent", "linkinvent")
COLORS = ((0.36, 0.68, 0.93), (0.38, 0.76, 0.50), (1.0, 0.70, 0.35))


@lru_cache(maxsize=16)
def ligand(target: str) -> Chem.Mol:
    if target not in atlas.CODES:
        raise ValueError(f"未知靶点：{target}")
    return atlas.load_ligand(target)


@lru_cache(maxsize=1)
def models() -> dict:
    return {
        mode: torch.load(
            atlas.ROOT / "REINVENT4" / "priors" / f"{mode}_transformer_pubchem.prior",
            map_location="cpu",
            weights_only=False,
        )
        for mode in MODES
    }


def eligible_bonds(mol: Chem.Mol) -> set[tuple[int, int]]:
    return {
        tuple(sorted((bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1)))
        for bond in mol.GetBonds()
        if not bond.IsInRing() and bond.GetBondType() == Chem.BondType.SINGLE
    }


def cuts_for(target: str, mode: str, raw_cuts: object, *, complete: bool) -> list[tuple[int, int]]:
    if mode not in MODES:
        raise ValueError(f"未知任务类型：{mode}")
    if not isinstance(raw_cuts, list):
        raise ValueError("切割键必须是列表")
    required = 1 if mode == "libinvent" else 2
    if len(raw_cuts) > required or (complete and len(raw_cuts) != required):
        raise ValueError(f"{mode} 需要选择 {required} 条切割键")
    allowed = eligible_bonds(ligand(target))
    cuts = []
    for raw in raw_cuts:
        if not isinstance(raw, list) or len(raw) != 2 or any(type(n) is not int for n in raw):
            raise ValueError("每条切割键须由两个原始 MOL2 原子 ID 表示")
        bond = tuple(sorted(raw))
        if bond not in allowed:
            raise ValueError(f"{raw} 不是可切割的非环单键")
        if bond in cuts:
            raise ValueError("同一条键不能重复选择")
        cuts.append(bond)
    return cuts


def fragments_for(target: str, cuts: list[tuple[int, int]]) -> list[Chem.Mol]:
    return atlas.split(ligand(target), cuts) if cuts else []


def fragment_info(fragment: Chem.Mol, *, include_svg: bool = True) -> dict:
    result = {
        "atom_ids": sorted(atlas.ids(fragment)),
        "heavy_atoms": fragment.GetNumHeavyAtoms(),
        "attachment_points": len(atlas.dummies(fragment)),
        "smiles": atlas.smi(fragment),
    }
    if include_svg:
        result["svg"] = draw_svg(fragment, width=410, height=260, annotate_ids=False)
    return result


def preview(target: str, mode: str, raw_cuts: object, raw_retained: object = None) -> dict:
    cuts = cuts_for(target, mode, raw_cuts, complete=False)
    required = 1 if mode == "libinvent" else 2
    fragments = fragments_for(target, cuts)
    expected = len(cuts) + 1 if cuts else 0
    if fragments and len(fragments) != expected:
        raise ValueError("切割后未得到预期数量的独立组分；请更换切割键")
    retained = [] if raw_retained is None else raw_retained
    if not isinstance(retained, list) or len(retained) > required:
        raise ValueError("保留组分数量无效")
    kept = _chosen_fragments(fragments, retained, len(retained))
    components = [atlas.ids(fragment) for fragment in kept]
    if len(cuts) == required and len(kept) == required:
        components.append(atlas.ids(next(fragment for fragment in fragments if fragment not in kept)))
    return {
        "cuts": [list(cut) for cut in cuts],
        "complete": len(cuts) == required,
        "fragments": [fragment_info(fragment) for fragment in fragments],
        "svg": draw_svg(ligand(target), cuts=tuple(cuts), components=tuple(components), interactive=True),
    }


def draw_svg(
    original: Chem.Mol,
    *,
    width: int = 1100,
    height: int = 680,
    cuts: tuple[tuple[int, int], ...] = (),
    components: tuple[set[int], ...] = (),
    interactive: bool = False,
    annotate_ids: bool = True,
) -> str:
    mol = Chem.Mol(original)
    if annotate_ids:
        for atom in mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIntProp("source_id")))
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.padding = 0.08
    options.minFontSize = 13
    options.maxFontSize = 28
    options.annotationFontScale = 0.62
    atom_colors = {}
    for index, component in enumerate(components):
        for atom_id in component:
            atom_colors[atom_id - 1] = COLORS[index]
    bond_colors = {}
    for a, b in cuts:
        bond = mol.GetBondBetweenAtoms(a - 1, b - 1)
        bond_colors[bond.GetIdx()] = (0.92, 0.19, 0.20)
    drawer.DrawMolecule(
        mol,
        highlightAtoms=list(atom_colors),
        highlightBonds=list(bond_colors),
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    if interactive:
        hit_lines = []
        for a, b in sorted(eligible_bonds(mol)):
            p = drawer.GetDrawCoords(a - 1)
            q = drawer.GetDrawCoords(b - 1)
            hit_lines.append(
                f'<line class="bond-hit" data-a="{a}" data-b="{b}" '
                f'x1="{p.x:.2f}" y1="{p.y:.2f}" x2="{q.x:.2f}" y2="{q.y:.2f}" '
                f'stroke="transparent" stroke-width="20" stroke-linecap="round" '
                f'pointer-events="stroke"><title>切割键 {a}—{b}</title></line>'
            )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("</svg>", "".join(hit_lines) + "</svg>")
        return svg[svg.index("<svg"):]
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg[svg.index("<svg"):]


def target_info(target: str) -> dict:
    mol = ligand(target)
    return {
        "target": target,
        "pdb": atlas.PDB[target],
        "ccd": atlas.CODES[target],
        "source_sha256": atlas.RAW_SHA256[target],
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "svg": draw_svg(mol, interactive=True),
    }


def _chosen_fragments(fragments: list[Chem.Mol], raw_retained: object, count: int) -> list[Chem.Mol]:
    if not isinstance(raw_retained, list) or len(raw_retained) != count:
        raise ValueError(f"请选择 {count} 个保留组分")
    selected = []
    for atom_ids in raw_retained:
        if not isinstance(atom_ids, list) or any(type(n) is not int for n in atom_ids):
            raise ValueError("保留组分须由原始 MOL2 原子 ID 列表表示")
        matched = [fragment for fragment in fragments if atlas.ids(fragment) == set(atom_ids)]
        if len(matched) != 1 or len(atom_ids) != len(set(atom_ids)):
            raise ValueError("所选组分与切分结果不一致")
        if matched[0] in selected:
            raise ValueError("不能重复选择同一组分")
        selected.append(matched[0])
    return selected


def _check_model_text(text: str, model: dict, label: str) -> int:
    try:
        return atlas.check_vocab(text, model, label)
    except AssertionError as error:
        raise ValueError(f"{label} 未通过 prior 词表或长度检查：{error}") from error


def _assembly_status(mol: Chem.Mol, joined: Chem.Mol | None) -> dict:
    if joined is None:
        raise ValueError("参考片段无法用 REINVENT4 的拼接逻辑复原")
    clean_joined = Chem.Mol(joined)
    for atom in clean_joined.GetAtoms():
        atom.SetAtomMapNum(0)
    connected = Chem.MolToSmiles(clean_joined, isomericSmiles=False) == Chem.MolToSmiles(mol, isomericSmiles=False)
    if not connected:
        raise ValueError("参考片段拼接后连接关系与原配体不一致")
    return {"connectivity_match": True, "stereochemistry_match": atlas.clean(joined) == atlas.clean(mol)}


def build_task(target: str, mode: str, spec: object, prior: dict) -> dict:
    if not isinstance(spec, dict):
        raise ValueError("任务设计必须是对象")
    cuts = cuts_for(target, mode, spec.get("cuts"), complete=True)
    mol = ligand(target)
    fragments = fragments_for(target, cuts)
    expected = 2 if mode == "libinvent" else 3
    if len(fragments) != expected:
        raise ValueError(f"切割后须得到 {expected} 个独立组分")
    kept = _chosen_fragments(fragments, spec.get("retained"), expected - 1)
    generated = next(fragment for fragment in fragments if fragment not in kept)
    if any(len(Chem.GetMolFrags(fragment)) != 1 for fragment in fragments):
        raise ValueError("切割后的组分不连续")
    if len(atlas.dummies(generated)) != expected - 1 or any(len(atlas.dummies(f)) != 1 for f in kept):
        raise ValueError("保留端和待生成区域的连接点数量不符合任务要求")
    if set.union(*(atlas.ids(f) for f in fragments)) != atlas.ids(mol):
        raise ValueError("切分未覆盖原配体全部原子")

    if mode == "libinvent":
        ordered_kept = kept
        input_smiles = atlas.conversions.convert_to_standardized_smiles(atlas.smi(kept[0]))
        output_smiles = atlas.smi(generated)
        joined = atlas.bond_maker.join_scaffolds_and_decorations(
            atlas.attachment_points.add_attachment_point_numbers(atlas.smi(kept[0]), canonicalize=False),
            output_smiles,
            keep_labels_on_atoms=True,
        )
    else:
        labeled_linker = Chem.Mol(generated)
        for atom in atlas.dummies(labeled_linker):
            atom.SetAtomMapNum(atom.GetIntProp("cut_id"))
        parsed = Chem.MolFromSmiles(atlas.smi(labeled_linker))
        cut_order = [atom.GetAtomMapNum() for atom in atlas.dummies(parsed)]
        ordered_kept = [
            next(f for f in kept if atlas.dummies(f)[0].GetIntProp("cut_id") == cut_id)
            for cut_id in cut_order
        ]
        for atom in parsed.GetAtoms():
            atom.SetAtomMapNum(0)
        output_smiles = Chem.MolToSmiles(parsed, canonical=False, isomericSmiles=True)
        input_smiles = "|".join(
            atlas.conversions.convert_to_standardized_smiles(atlas.smi(fragment))
            for fragment in ordered_kept
        )
        joined = atlas.bond_maker.join_scaffolds_and_decorations(
            atlas.attachment_points.add_attachment_point_numbers(output_smiles, canonicalize=False),
            "|".join(atlas.smi(fragment) for fragment in ordered_kept),
            keep_labels_on_atoms=True,
        )

    assembly = _assembly_status(mol, joined)
    input_parts = input_smiles.split("|")
    for fragment, standardized in zip(ordered_kept, input_parts, strict=True):
        parsed_input = Chem.MolFromSmiles(standardized)
        if parsed_input is None or atlas.clean(parsed_input) != atlas.clean(fragment):
            raise ValueError("模型输入标准化改变了保留组分")
    input_tokens = _check_model_text(input_smiles, prior, f"{target} {mode} 输入")
    output_tokens = _check_model_text(output_smiles, prior, f"{target} {mode} 参考输出")
    note = spec.get("note", "")
    if not isinstance(note, str) or len(note) > 4000:
        raise ValueError("备注须为不超过 4000 字的文本")
    return {
        "status": "designed",
        "cuts": [list(cut) for cut in cuts],
        "retained": [fragment_info(f, include_svg=False) for f in ordered_kept],
        "generated_reference": fragment_info(generated, include_svg=False),
        "input_smiles": input_smiles,
        "reference_output_smiles": output_smiles,
        "input_tokens": input_tokens,
        "reference_output_tokens": output_tokens,
        "validation": assembly,
        "note": note,
    }


def task_record(target: str, mode: str, spec: object) -> dict:
    if not isinstance(spec, dict):
        raise ValueError(f"{target} {mode} 的任务设计必须是对象")
    note = spec.get("note", "")
    if not isinstance(note, str) or len(note) > 4000:
        raise ValueError(f"{target} {mode} 的备注须为不超过 4000 字的文本")
    skipped = spec.get("skip", False)
    if type(skipped) is not bool:
        raise ValueError(f"{target} {mode} 的跳过标记必须是布尔值")
    cuts = cuts_for(target, mode, spec.get("cuts", []), complete=False)
    raw_retained = spec.get("retained", [])
    required = 1 if mode == "libinvent" else 2
    if not isinstance(raw_retained, list) or len(raw_retained) > required:
        raise ValueError(f"{target} {mode} 的保留组分数量无效")
    fragments = fragments_for(target, cuts)
    if fragments and len(fragments) != len(cuts) + 1:
        raise ValueError(f"{target} {mode} 切割后未得到预期数量的组分")
    retained = _chosen_fragments(fragments, raw_retained, len(raw_retained))
    if any(len(atlas.dummies(fragment)) != 1 for fragment in retained):
        raise ValueError(f"{target} {mode} 的保留组分必须各有一个连接点")
    if not skipped and len(cuts) == required and len(retained) == required:
        return build_task(target, mode, spec, models()[mode])
    status = "skipped" if skipped else "incomplete" if cuts or retained or note else "not_started"
    return {
        "status": status,
        "cuts": [list(cut) for cut in cuts],
        "retained": [sorted(atlas.ids(fragment)) for fragment in retained],
        "note": note,
    }


def build_snapshot(raw: object) -> dict:
    if not isinstance(raw, dict) or set(raw) - set(atlas.CODES):
        raise ValueError("设计数据必须是对象，且只能包含本项目的 16 个靶点")
    targets = {}
    for target in atlas.CODES:
        source = atlas.ROOT / "real-world_dataset" / target / "crystal.mol2"
        if hashlib.sha256(source.read_bytes()).hexdigest() != atlas.RAW_SHA256[target]:
            raise ValueError(f"{target} 原始 MOL2 已改变；请重启界面并重新核对设计")
        tasks = raw.get(target, {})
        if not isinstance(tasks, dict) or set(tasks) - set(MODES):
            raise ValueError(f"{target} 的任务数据必须只包含 LibINVENT/LinkINVENT")
        targets[target] = {
            "pdb": atlas.PDB[target],
            "ccd": atlas.CODES[target],
            "source_mol2": f"real-world_dataset/{target}/crystal.mol2",
            "source_sha256": atlas.RAW_SHA256[target],
            "reference_ligand_smiles": atlas.smi(ligand(target)),
            "tasks": {
                mode: task_record(target, mode, tasks.get(mode, {})) for mode in MODES
            },
        }
    now = datetime.now(ZoneInfo("Asia/Shanghai"))
    counts = {
        status: sum(task["status"] == status for info in targets.values() for task in info["tasks"].values())
        for status in ("designed", "skipped", "incomplete", "not_started")
    }
    return {
        "schema_version": 2,
        "created_at": now.isoformat(timespec="microseconds"),
        "source": "manual bond and retained-component selection",
        "rdkit_version": rdBase.rdkitVersion,
        "model_ids": ({mode: models()[mode]["metadata"]["model_id"] for mode in MODES}
                      if counts["designed"] else {}),
        "status_counts": counts,
        "targets": targets,
    }


def report_html(snapshot: dict) -> str:
    rows = []
    for target, info in snapshot["targets"].items():
        tasks = []
        for mode, label in (("libinvent", "LibINVENT"), ("linkinvent", "LinkINVENT")):
            task = info["tasks"][mode]
            cuts = ", ".join(f"{a}—{b}" for a, b in task["cuts"])
            status_label = {"designed": "已设计", "skipped": "跳过", "incomplete": "未完成", "not_started": "尚未开始"}[task["status"]]
            if task["status"] == "designed":
                retained = "; ".join(", ".join(map(str, frag["atom_ids"])) for frag in task["retained"])
                stereo = "一致" if task["validation"]["stereochemistry_match"] else "不一致，需审阅"
                details = (
                    f'<dt>保留原子 ID</dt><dd>{retained}</dd>'
                    f'<dt>待生成原子 ID</dt><dd>{", ".join(map(str, task["generated_reference"]["atom_ids"]))}</dd>'
                    f'<dt>模型输入 SMILES</dt><dd><code>{escape(task["input_smiles"])}</code></dd>'
                    f'<dt>参考输出 SMILES</dt><dd><code>{escape(task["reference_output_smiles"])}</code></dd>'
                    f'<dt>拼接立体化学</dt><dd>{stereo}</dd>'
                )
            else:
                retained = "; ".join(", ".join(map(str, ids)) for ids in task["retained"])
                details = f'<dt>已选保留原子 ID</dt><dd>{retained or "—"}</dd>'
            tasks.append(
                f'<section><h3>{label} · {status_label}</h3><img src="figures/{target}_{mode}.svg" alt="{target} {label} 切分图">'
                f'<dl><dt>切割键</dt><dd>{cuts or "—"}</dd>{details}'
                f'<dt>备注</dt><dd>{escape(task["note"]) or "—"}</dd></dl></section>'
            )
        rows.append(
            f'<article id="{target}"><h2>{target.upper()} · {info["pdb"]} / {info["ccd"]}</h2>'
            f'<p>原始 MOL2 SHA-256：<code>{info["source_sha256"]}</code></p>'
            + "".join(tasks) + "</article>"
        )
    links = " ".join(f'<a href="#{target}">{target.upper()}</a>' for target in snapshot["targets"])
    return (
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>手动任务设计记录</title><style>'
        'body{font:16px/1.6 system-ui,sans-serif;color:#183047;background:#edf3f7;margin:0}'
        'header,main{max-width:1250px;margin:auto;padding:24px}header,article{background:white}'
        'article{margin:20px 0;padding:24px;border-radius:12px}section{border-top:1px solid #d9e3ea;padding:15px 0}'
        'img{display:block;max-width:100%;height:auto;background:#fff;border:1px solid #d9e3ea}'
        'dl{display:grid;grid-template-columns:180px 1fr;gap:7px 14px}dt{font-weight:700}'
        'dd{margin:0;overflow-wrap:anywhere}code{font-size:13px}nav{display:flex;flex-wrap:wrap;gap:12px}'
        'a{color:#09639c}@media(max-width:650px){dl{display:block}dt{margin-top:12px}}'
        '</style><header><h1>16 个靶点的手动任务设计</h1>'
        f'<p>记录时间：{escape(snapshot["created_at"])}；已设计 {snapshot["status_counts"]["designed"]} 项，'
        f'跳过 {snapshot["status_counts"]["skipped"]} 项。完整数据见 <a href="designs.json">designs.json</a>。'
        '已设计任务中橙色为参考配体的待生成区域；红色为切割键。图上的数字是原始 MOL2 原子 ID。</p>'
        f'<nav>{links}</nav></header><main>{"".join(rows)}</main></html>'
    )


def save_snapshot(raw: object, output_root: Path = HERE) -> Path:
    snapshot = build_snapshot(raw)
    figures = {}
    for target, info in snapshot["targets"].items():
        for mode in MODES:
            task = info["tasks"][mode]
            if task["status"] == "designed":
                components = tuple([set(fragment["atom_ids"]) for fragment in task["retained"]]
                                   + [set(task["generated_reference"]["atom_ids"])])
            else:
                components = tuple(set(ids) for ids in task["retained"])
            figures[f"{target}_{mode}.svg"] = draw_svg(
                ligand(target),
                cuts=tuple(tuple(cut) for cut in task["cuts"]),
                components=components,
            )
    timestamp = datetime.fromisoformat(snapshot["created_at"]).strftime("%Y-%m-%dT%H-%M-%S-%f%z")
    folder = output_root / timestamp
    folder.mkdir(parents=True, exist_ok=False)
    (folder / "figures").mkdir()
    (folder / "designs.json").write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (folder / "index.html").write_text(report_html(snapshot), encoding="utf-8")
    for name, svg in figures.items():
        (folder / "figures" / name).write_text(svg, encoding="utf-8")
    return folder
