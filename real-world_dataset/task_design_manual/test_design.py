"""Regression checks for manual design validation and immutable records."""

from pathlib import Path
import json
import sys
import xml.etree.ElementTree as ET

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import design
import app


def atlas_examples() -> dict:
    selections = {}
    for target in design.atlas.CODES:
        lib_cut, seed, link_cuts, _, _ = design.atlas.DESIGNS[target]
        lib_parts = design.fragments_for(target, [tuple(lib_cut)])
        core = next(part for part in lib_parts if seed in design.atlas.ids(part))
        link_parts = design.fragments_for(target, [tuple(cut) for cut in link_cuts])
        ends = [part for part in link_parts if len(design.atlas.dummies(part)) == 1]
        selections[target] = {
            "libinvent": {"cuts": [list(lib_cut)], "retained": [sorted(design.atlas.ids(core))], "note": ""},
            "linkinvent": {"cuts": [list(cut) for cut in link_cuts],
                           "retained": [sorted(design.atlas.ids(end)) for end in ends], "note": ""},
        }
    return selections


def test_all_atlas_examples_and_distinct_records(tmp_path):
    selections = atlas_examples()
    snapshot = design.build_snapshot(selections)
    assert len(snapshot["targets"]) == 16
    assert all(len(item["tasks"]) == 2 for item in snapshot["targets"].values())
    assert snapshot["schema_version"] == 2
    assert snapshot["status_counts"]["designed"] == 32
    assert snapshot["targets"]["dyr"]["tasks"]["linkinvent"]["validation"]["stereochemistry_match"] is False
    assert all(task["validation"]["connectivity_match"]
               for item in snapshot["targets"].values() for task in item["tasks"].values())

    first = design.save_snapshot(selections, tmp_path)
    second = design.save_snapshot(selections, tmp_path)
    assert first != second
    assert (first / "designs.json").is_file()
    assert (first / "index.html").is_file()
    assert len(list((first / "figures").glob("*.svg"))) == 32
    assert len(list((second / "figures").glob("*.svg"))) == 32
    assert json.loads((first / "designs.json").read_text())["targets"]["aa2ar"]["tasks"]["libinvent"]["cuts"] == [[9, 10]]
    assert (first / "index.html").read_text().count('src="figures/') == 32
    ET.fromstring((first / "figures" / "aa2ar_libinvent.svg").read_text())
    assert app.RECORD_PATH.fullmatch(f"/{first.name}/figures/aa2ar_libinvent.svg")


def test_invalid_cut_and_component_are_rejected():
    selections = atlas_examples()
    with pytest.raises(ValueError, match="非环单键"):
        design.preview("aa2ar", "libinvent", [[1, 2]])
    selections["aa2ar"]["libinvent"]["retained"] = [[1]]
    with pytest.raises(ValueError, match="不一致"):
        design.build_task("aa2ar", "libinvent", selections["aa2ar"]["libinvent"], design.models()["libinvent"])


def test_partial_record_preserves_skipped_incomplete_and_empty_tasks(tmp_path):
    selections = atlas_examples()
    partial = {
        "aa2ar": {
            "libinvent": selections["aa2ar"]["libinvent"],
            "linkinvent": {"skip": True, "note": "没有合适的两个保留端"},
        },
        "abl1": {"libinvent": {"cuts": selections["abl1"]["libinvent"]["cuts"], "note": "待选择保留端"}},
    }
    snapshot = design.build_snapshot(partial)
    assert snapshot["status_counts"] == {
        "designed": 1, "skipped": 1, "incomplete": 1, "not_started": 29,
    }
    assert snapshot["targets"]["aa2ar"]["tasks"]["linkinvent"]["note"] == "没有合适的两个保留端"
    assert snapshot["targets"]["abl1"]["tasks"]["libinvent"]["cuts"]
    assert snapshot["targets"]["ppara"]["tasks"]["linkinvent"]["status"] == "not_started"
    folder = design.save_snapshot(partial, tmp_path)
    assert len(list((folder / "figures").glob("*.svg"))) == 32
    html = (folder / "index.html").read_text()
    assert "跳过" in html and "未完成" in html and "尚未开始" in html


def test_empty_record_is_allowed_and_bad_finished_task_is_rejected(tmp_path):
    assert design.build_snapshot({})["status_counts"]["not_started"] == 32
    selections = atlas_examples()
    selections["aa2ar"]["libinvent"]["retained"] = [[1]]
    with pytest.raises(ValueError, match="不一致"):
        design.save_snapshot({"aa2ar": {"libinvent": selections["aa2ar"]["libinvent"]}}, tmp_path)
    assert not list(tmp_path.iterdir())


def test_interactive_svg_contains_clickable_bonds():
    svg = design.target_info("aa2ar")["svg"]
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert 'data-a="9" data-b="10"' in svg
