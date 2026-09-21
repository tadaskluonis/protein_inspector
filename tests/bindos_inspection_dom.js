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
    annotationRows[a.annotation_id] = {hidden: false, style: {}, dataset: {annotation: a.annotation_id},
        addEventListener() {}};
}
const detailsNode = {textContent: ''};

function bySelector(sel) {
    let m = sel.match(/data-layer="([^"]+)"/);
    if (m) return layerBoxes[m[1]] || null;
    m = sel.match(/data-annotation="([^"]+)"/);
    if (m) return annotationRows[m[1]] || null;
    return null;
}

const rafQueue = [];
let frame = 0;
global.window = {};
global.CSS = {escape: (s) => s};
global.requestAnimationFrame = (fn) => rafQueue.push(fn);
global.document = {
    getElementById: (id) => (id === 'bindos-inspection-state'
        ? {textContent: stateMatch[1]}
        : (id === 'bindos-residue-details' ? detailsNode : null)),
    querySelector: bySelector,
    querySelectorAll: (sel) => (sel === '[data-annotation]'
        ? Object.values(annotationRows)
        : (sel === '[data-layer]' ? Object.values(layerBoxes) : [])),
    addEventListener: () => {},
};

// --- stub renderer, published the way src/parts/embed.js publishes it -----
const viewerId = (state.viewer && state.viewer.config && state.viewer.config.viewer_id) || 'stub';
const object = {};
let renders = 0;
const renderer = {
    residueNumbers, chains,
    currentObjectName: 'target',
    objectsData: {target: object},
    colorsNeedUpdate: false, plddtColorsNeedUpdate: false,
    residueSelection: new Set(),
    setResidueSelection(next) { this.residueSelection = next; },
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
for (let i = 0; i < 700 && rafQueue.length; i++) {
    if (++frame === rafDelay) publishViewer();
    rafQueue.splice(0, rafQueue.length).forEach((fn) => fn());
}

const paint = () => (object.color && object.color.type === 'advanced'
    ? object.color.value.position : null);
const report = {initial: paint(), renders, frames: frame, steps: {}};

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
report.finalRenders = renders;
process.stdout.write(JSON.stringify(report));
