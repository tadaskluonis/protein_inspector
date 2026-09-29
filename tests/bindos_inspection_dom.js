// Execute the inspection bundle's behaviour script the way a browser would,
// against a stub DOM and a stub renderer, and report what it did.
//
//   node tests/bindos_inspection_dom.js <inspection.html> '<spec json>'
//
// spec: {"chains": ["A", ...], "residueNumbers": [1, ...], "rafDelay": 0}
//
// WHY THIS EXISTS. The Python tests assert that strings such as
// "syncVisibleLayers" appear in the emitted HTML. That is satisfied by a
// handler that never runs, and for a while one did not: the script looked the
// renderer up through a bare `viewerApi`, which is `let viewerApi` inside
// src/app/main.js and therefore not reachable from a separate <script>. The
// page's real contract is window.py2dmol_viewers[id] = {renderer}, set in
// src/parts/embed.js. SO THIS HARNESS DELIBERATELY DOES NOT DEFINE viewerApi.
// Defining it would re-hide the exact bug this file was written to catch.
'use strict';
const fs = require('fs');

const html = fs.readFileSync(process.argv[2], 'utf8');
const spec = JSON.parse(process.argv[3] || '{}');
const chains = spec.chains || [];
const residueNumbers = spec.residueNumbers || [];
const rafDelay = spec.rafDelay || 0;

// --- pull the two scripts out of the page --------------------------------
const stateMatch = html.match(
    /<script id="bindos-inspection-state" type="application\/json">([\s\S]*?)<\/script>/);
if (!stateMatch) throw new Error('no bindos-inspection-state in page');
const state = JSON.parse(stateMatch[1].replace(/<\\\//g, '</'));

// py2Dmol's own loader shim also opens with `<script>(function(){`, so match
// on the handler name rather than on the shape of the opening bytes.
const blocks = html.match(/<script>[\s\S]*?<\/script>/g) || [];
const behaviour = blocks.filter((b) => b.includes('syncVisibleLayers'));
if (behaviour.length !== 1) {
    throw new Error(`expected exactly 1 behaviour script, found ${behaviour.length}`);
}
const source = behaviour[0].replace(/^<script>/, '').replace(/<\/script>$/, '');

// --- stub DOM -------------------------------------------------------------
const layerBoxes = {};
for (const layer of state.layers) {
    layerBoxes[layer.layer_id] = {checked: true, handlers: [], hidden: false,
        addEventListener(type, fn) { this.handlers.push(fn); }};
}
const annotationRows = {};
for (const a of state.annotations) {
    annotationRows[a.annotation_id] = {hidden: false, style: {}, attrs: {},
        dataset: {annotation: a.annotation_id},
        addEventListener(type, fn) { (this.handlers = this.handlers || []).push(fn); },
        setAttribute(k, v) { this.attrs[k] = v; },
        removeAttribute(k) { delete this.attrs[k]; },
        click() { (this.handlers || []).forEach((fn) => fn()); }};
}
const detailsNode = {textContent: ''};

// The panel's three buttons. They are looked up with getElementById and the
// script guards each one, so a harness that omitted them would let a dead
// handler pass -- the same failure mode this file exists to catch.
function stubButton() {
    return {handlers: [], attrs: {}, textContent: '',
            addEventListener(type, fn) { this.handlers.push(fn); },
            setAttribute(k, v) { this.attrs[k] = v; },
            getAttribute(k) { return this.attrs[k]; },
            removeAttribute(k) { delete this.attrs[k]; },
            click() { this.handlers.forEach((fn) => fn()); }};
}
const buttons = {
    'bindos-layers-all': stubButton(),
    'bindos-layers-none': stubButton(),
    'bindos-clear-residue': stubButton(),
    'bindos-controls': stubButton(),
    'bindos-morph-note': {textContent: ''},
};
// The stage carries the controls-collapsed state as an attribute; without it
// the toggle silently does nothing and the harness would not notice.
const stage = {attrs: {'data-controls': '1'}, clientWidth: spec.stageWidth || 1200,
    setAttribute(k, v) { this.attrs[k] = v; },
    getAttribute(k) { return this.attrs[k]; },
    querySelector() { return null; }};
buttons['bindos-stage'] = stage;

// One button per conformation, discovered the way the page discovers them.
const morphButtons = ((state.morph && state.morph.conformers) || []).map((c, i) => {
    const b = stubButton();
    b.attrs['data-conf'] = String(i);
    return b;
});

function bySelector(sel) {
    let m = sel.match(/data-layer="([^"]+)"/);
    if (m) return layerBoxes[m[1]] || null;
    m = sel.match(/data-annotation="([^"]+)"/);
    if (m) return annotationRows[m[1]] || null;
    return null;
}

const rafQueue = [];
let frame = 0;
// The page listens on window (resize -> control autofit) and schedules work
// with requestAnimationFrame; a bare {} makes the whole behaviour script die
// on the first addEventListener, which is a stub gap and not a page bug.
// A clock the drain loop drives, not the wall clock: the animation is written
// against elapsed milliseconds, and a harness that drains 700 frames in 3 ms
// would never reach the end of one.
let clock = 0;
global.window = {addEventListener() {}, removeEventListener() {},
    performance: {now: () => clock}};
global.requestAnimationFrame = global.requestAnimationFrame || (fn => setTimeout(fn, 0));
global.CSS = {escape: (s) => s};
global.requestAnimationFrame = (fn) => rafQueue.push(fn);
global.document = {
    getElementById: (id) => (id === 'bindos-inspection-state'
        ? {textContent: stateMatch[1]}
        : (id === 'bindos-residue-details' ? detailsNode : (buttons[id] || null))),
    querySelector: bySelector,
    querySelectorAll: (sel) => (sel === '[data-annotation]'
        ? Object.values(annotationRows)
        : (sel === '[data-layer]' ? Object.values(layerBoxes)
            : (sel === '.bm-btn' ? morphButtons : []))),
    addEventListener: () => {},
};

// --- stub renderer, published the way src/parts/embed.js publishes it -----
const viewerId = (state.viewer && state.viewer.config && state.viewer.config.viewer_id) || 'stub';
const object = {};
// ENDPOINTS ONLY, as the file ships them: one frame per conformation, each a
// distinct coordinate set so an interpolated frame is distinguishable from a
// copied one. The page is expected to expand these in place.
if (state.morph && state.morph.mode === 'browser') {
    const n = Math.max(1, residueNumbers.length);
    object.frames = state.morph.conformers.map((c, k) => ({
        name: c.label,
        coords: Array.from({length: n}, (_, i) => [i, k * 10, 0]),
    }));
}
let renders = 0;
let uiUpdates = 0;
const renderer = {
    residueNumbers, chains,
    currentObjectName: 'target',
    currentFrame: 0,
    objectsData: {target: object},
    colorsNeedUpdate: false, plddtColorsNeedUpdate: false,
    residueSelection: new Set(),
    framesVisited: [],
    setResidueSelection(next) { this.residueSelection = next; },
    setFrame(i) { this.currentFrame = i; this.framesVisited.push(i); renders++; },
    updateUIControls() { uiUpdates++; },
    render() { renders++; },
};
function publishViewer() {
    window.py2dmol_viewers = {[viewerId]: {renderer}};
    window.viewerConfig = {viewer_id: viewerId};
}
if (rafDelay === 0) publishViewer();

// --- run it ---------------------------------------------------------------
new Function(source)();

// Drain animation frames; the viewer appears late when rafDelay > 0, which is
// what happens in the browser (the bundle builds the renderer on load).
function drain(limit) {
    for (let i = 0; i < limit && rafQueue.length; i++) {
        if (++frame === rafDelay) publishViewer();
        clock += 16;
        rafQueue.splice(0, rafQueue.length).forEach((fn) => fn());
    }
}
drain(700);

const paint = () => (object.color && object.color.type === 'advanced'
    ? object.color.value.position : null);
const report = {initial: paint(), renders, frames: frame, steps: {}};

// --- morph ----------------------------------------------------------------
// The file ships N endpoints; the page must expand them to N + (N-1)*(steps-1)
// + ... frames on load and animate THROUGH them on a button press, landing on
// the target conformation's own frame.
if (state.morph && state.morph.mode === 'browser' && morphButtons.length) {
    const target = morphButtons.length - 1;
    const stride = state.morph.steps_between;
    report.morph = {
        framesInFile: state.morph.conformers.length,
        framesAfterExpansion: object.frames.length,
        uiUpdates,
        pressedBefore: morphButtons.map((b) => b.getAttribute('aria-pressed')),
    };
    renderer.framesVisited.length = 0;
    morphButtons[target].click();
    drain(400);
    const visited = renderer.framesVisited;
    report.morph.visitedCount = visited.length;
    report.morph.firstFrame = visited[0];
    report.morph.landedOn = visited[visited.length - 1];
    report.morph.expectedLanding = target * stride;
    report.morph.monotonic = visited.every((v, i) => i === 0 || v >= visited[i - 1]);
    report.morph.pressedAfter = morphButtons.map((b) => b.getAttribute('aria-pressed'));
    report.morph.note = buttons['bindos-morph-note'].textContent;
    // An interpolated frame must be a genuine blend, not a copy of an endpoint.
    const mid = object.frames[Math.floor(stride / 2)];
    report.morph.midpointY = mid ? mid.coords[0][1] : null;
}

// --- controls collapse ----------------------------------------------------
buttons['bindos-controls'].click();
report.controlsAfterClick = stage.getAttribute('data-controls');

for (const step of (spec.steps || [])) {
    for (const [layerId, checked] of Object.entries(step.set)) {
        layerBoxes[layerId].checked = checked;
    }
    const box = layerBoxes[Object.keys(step.set)[0]];
    box.handlers.forEach((fn) => fn());
    report.steps[step.name] = {
        paint: paint(),
        hiddenRows: Object.entries(annotationRows).filter(([, r]) => r.hidden).length,
        selection: renderer.residueSelection.size,
    };
}
// --- panel buttons --------------------------------------------------------
// "Untick all" must clear every checkbox and repaint to no colour; "Tick all"
// must restore it. "Clear" must empty the detail pane and drop the selection.
if (spec.buttons !== false) {
    buttons['bindos-layers-none'].click();
    report.untickAll = {
        paint: paint(),
        checked: Object.values(layerBoxes).filter((b) => b.checked).length,
        hiddenRows: Object.entries(annotationRows).filter(([, r]) => r.hidden).length,
    };
    buttons['bindos-layers-all'].click();
    report.tickAll = {
        paint: paint(),
        checked: Object.values(layerBoxes).filter((b) => b.checked).length,
    };
    const someAnnotation = state.annotations.find((a) => a.resolved);
    if (someAnnotation) window.bindosInspection.selectAnnotation(someAnnotation.annotation_id);
    report.afterSelect = {details: detailsNode.textContent.length,
                          selection: renderer.residueSelection.size};
    buttons['bindos-clear-residue'].click();
    report.afterClear = {details: detailsNode.textContent,
                         selection: renderer.residueSelection.size};
}

report.finalRenders = renders;
process.stdout.write(JSON.stringify(report));
