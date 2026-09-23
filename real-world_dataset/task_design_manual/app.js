"use strict";

const DRAFT_KEY = "task-design-manual-v1";
const MODES = ["libinvent", "linkinvent"];
let targets = [];
let currentTarget = null;
let currentMode = "libinvent";
let currentInfo = null;
let designs = {};
let previewNumber = 0;

const byId = id => document.getElementById(id);
const task = () => designs[currentTarget][currentMode];
const needed = () => currentMode === "libinvent" ? 1 : 2;
const keyOf = ids => ids.join(",");

async function api(path, body) {
  const response = await fetch(path, body === undefined ? undefined : {
    method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify(body)
  });
  const result = await response.json();
  if (!response.ok) throw new Error(result.error || `请求失败：${response.status}`);
  return result;
}

function blankTask() { return {cuts: [], retained: [], note: ""}; }
function blankDesigns() {
  return Object.fromEntries(targets.map(item => [item.target, {
    libinvent: blankTask(), linkinvent: blankTask()
  }]));
}
function saveDraft() {
  try {
    localStorage.setItem(DRAFT_KEY, JSON.stringify({
      hashes: Object.fromEntries(targets.map(item => [item.target, item.sha256])), designs
    }));
  } catch (_) { /* Editing and recording still work when browser storage is disabled. */ }
  renderProgress();
}
function restoreDraft() {
  try {
    const saved = localStorage.getItem(DRAFT_KEY);
    if (!saved) return;
    const parsed = JSON.parse(saved);
    if (targets.every(item => parsed.hashes[item.target] === item.sha256) &&
        targets.every(item => MODES.every(mode => parsed.designs[item.target][mode]))) {
      designs = parsed.designs;
    }
  } catch (_) { /* Invalid or old browser draft: start with blank selections. */ }
}
function completed(item, mode) {
  const entry = designs[item][mode];
  const count = mode === "libinvent" ? 1 : 2;
  return entry.cuts.length === count && entry.retained.length === count;
}
function renderProgress() {
  const count = targets.reduce((sum, item) => sum + MODES.filter(mode => completed(item.target, mode)).length, 0);
  byId("progress").textContent = `${count} / 32 项完成`;
  byId("targets").replaceChildren(...targets.map(item => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = item.target === currentTarget ? "active" : "";
    button.innerHTML = `<span>${item.target.toUpperCase()}</span><span class="marks"><span class="${completed(item.target, "libinvent") ? "done" : ""}">●</span><span class="${completed(item.target, "linkinvent") ? "done" : ""}">●</span></span>`;
    button.addEventListener("click", () => showTarget(item.target));
    return button;
  }));
}
function message(text, kind = "info", html = false) {
  const panel = byId("message");
  panel.className = kind;
  if (html) panel.innerHTML = text; else panel.textContent = text;
}
function clearMessage() { byId("message").className = ""; byId("message").textContent = ""; }

async function showTarget(target) {
  currentTarget = target;
  currentInfo = null;
  renderProgress();
  byId("target-title").textContent = `${target.toUpperCase()} · 载入中…`;
  try {
    const info = await api(`/api/targets/${target}`);
    if (currentTarget !== target) return;
    currentInfo = info;
    byId("target-title").textContent = `${target.toUpperCase()} · ${info.pdb} / ${info.ccd}`;
    byId("source").textContent = `${info.heavy_atoms} 个重原子 · MOL2 SHA-256 ${info.source_sha256.slice(0, 14)}…`;
    byId("molecule").innerHTML = info.svg;
    byId("molecule").onclick = onBondClick;
    await showMode();
  } catch (error) { message(error.message, "error"); }
}

async function showMode() {
  if (!currentInfo) return;
  document.querySelectorAll(".tabs button").forEach(button => {
    button.classList.toggle("active", button.dataset.mode === currentMode);
  });
  byId("task-hint").textContent = currentMode === "libinvent"
    ? "选择一条切割键，再点击要固定保留的一个组分；另一组分为参考待生成区域。"
    : "选择两条切割键，再点击要固定保留的两个端组分；剩余组分为参考 linker。保留端 A/B 会按 linker 出口顺序记录。";
  byId("note").value = task().note;
  byId("cut-count").textContent = `${task().cuts.length} / ${needed()}`;
  byId("keep-count").textContent = `${task().retained.length} / ${needed()}`;
  renderCuts();
  await updatePreview();
}

function renderCuts() {
  const selected = new Set(task().cuts.map(cut => [...cut].sort((a, b) => a - b).join("-")));
  byId("molecule").querySelectorAll(".bond-hit").forEach(line => {
    line.classList.toggle("selected", selected.has(`${line.dataset.a}-${line.dataset.b}`));
  });
  byId("cut-list").textContent = task().cuts.length
    ? `已选切割键：${task().cuts.map(cut => cut.join("—")).join("；")}` : "尚未选择切割键";
  byId("cut-count").textContent = `${task().cuts.length} / ${needed()}`;
}

function onBondClick(event) {
  const line = event.target.closest(".bond-hit");
  if (!line || !currentInfo) return;
  const cut = [Number(line.dataset.a), Number(line.dataset.b)];
  const list = task().cuts;
  const index = list.findIndex(item => item[0] === cut[0] && item[1] === cut[1]);
  if (index >= 0) list.splice(index, 1);
  else if (list.length < needed()) list.push(cut);
  else { message(`请先取消一条切割键；${currentMode} 最多选择 ${needed()} 条。`, "error"); return; }
  task().retained = [];
  saveDraft();
  renderCuts();
  updatePreview();
}

async function updatePreview() {
  const serial = ++previewNumber;
  const entry = task();
  byId("components").replaceChildren();
  byId("keep-count").textContent = `${entry.retained.length} / ${needed()}`;
  if (!entry.cuts.length) {
    byId("molecule").innerHTML = currentInfo.svg;
    renderCuts();
    return;
  }
  try {
    const result = await api("/api/preview", {
      target: currentTarget, mode: currentMode, cuts: entry.cuts, retained: entry.retained
    });
    if (serial !== previewNumber) return;
    byId("molecule").innerHTML = result.svg;
    renderCuts();
    if (!result.complete) message("还需选择一条切割键，之后才能选择保留组分。", "info");
    else clearMessage();
    result.fragments.forEach((fragment, index) => {
      const card = document.createElement("div");
      const kept = entry.retained.some(ids => keyOf(ids) === keyOf(fragment.atom_ids));
      const generated = result.complete && entry.retained.length === needed() && !kept;
      card.className = `component${kept ? " retained" : ""}${generated ? " generated" : ""}`;
      card.innerHTML = `<div class="meta"><span>组分 ${index + 1} · ${fragment.heavy_atoms} 个重原子</span><span>${kept ? "保留" : generated ? "待生成" : `${fragment.attachment_points} 个连接点`}</span></div><div class="fragment-svg"></div><div class="smiles"></div>`;
      card.querySelector(".fragment-svg").innerHTML = fragment.svg;
      card.querySelector(".smiles").textContent = fragment.smiles;
      card.title = "原子 ID：" + fragment.atom_ids.join(", ");
      card.addEventListener("click", () => chooseFragment(fragment, result));
      byId("components").append(card);
    });
  } catch (error) { if (serial === previewNumber) message(error.message, "error"); }
}

function chooseFragment(fragment, result) {
  if (!result.complete) return;
  if (fragment.attachment_points !== 1) {
    message("保留端必须只有一个连接点；请选择两侧的端组分。", "error");
    return;
  }
  const list = task().retained;
  const index = list.findIndex(ids => keyOf(ids) === keyOf(fragment.atom_ids));
  if (index >= 0) list.splice(index, 1);
  else if (list.length < needed()) list.push(fragment.atom_ids);
  else { message(`请先取消一个保留组分；当前任务需要 ${needed()} 个。`, "error"); return; }
  saveDraft();
  updatePreview();
}

async function recordAll() {
  const missing = targets.flatMap(item => MODES.filter(mode => !completed(item.target, mode)).map(mode => `${item.target.toUpperCase()} ${mode}`));
  if (missing.length) {
    message(`仍有 ${missing.length} 项未完成：${missing.join("、")}`, "error");
    return;
  }
  const button = byId("record");
  button.disabled = true;
  message("正在核对 32 项设计并生成记录，请稍候…");
  try {
    const result = await api("/api/record", {designs});
    message(`记录已保存至 ${result.directory}。<a href="${result.report}" target="_blank" rel="noopener">打开报告</a>`, "success", true);
  } catch (error) { message(error.message, "error"); }
  finally { button.disabled = false; }
}

async function start() {
  try {
    targets = await api("/api/targets");
    designs = blankDesigns();
    restoreDraft();
    document.querySelectorAll(".tabs button").forEach(button => button.addEventListener("click", () => {
      currentMode = button.dataset.mode;
      showMode();
    }));
    byId("note").addEventListener("input", event => { task().note = event.target.value; saveDraft(); });
    byId("record").addEventListener("click", recordAll);
    await showTarget(targets[0].target);
  } catch (error) { message(error.message, "error"); }
}

start();
