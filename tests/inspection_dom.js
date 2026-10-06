// Execute the inspection bundle's behaviour script the way a browser would,
// against a stub DOM and a stub renderer, and report what it did.
//
//   node tests/inspection_dom.js <inspection.html> '<spec json>'
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
    /<script id="pinsp-inspection-state" type="application\/json">([\s\S]*?)<\/script>/);
if (!stateMatch) throw new Error('no pinsp-inspection-state in page');
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
// The residue card. It was absent from the harness, so every assertion about
// it was an assertion about `if (!card) return` -- and the card is now the
// answer to a click on an UNANNOTATED residue, which is the whole of what
// clicking the structure gets you on most of a structure.
const cardNode = {innerHTML: '', hidden: true};

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
    'pinsp-layers-all': stubButton(),
    'pinsp-layers-none': stubButton(),
    'pinsp-clear-residue': stubButton(),
    'pinsp-controls': stubButton(),
    'pinsp-morph-note': {textContent: ''},
};
// The stage carries the controls-collapsed state as an attribute; without it
// the toggle silently does nothing and the harness would not notice.
const stage = {attrs: {'data-controls': '1'}, clientWidth: spec.stageWidth || 1200,
    setAttribute(k, v) { this.attrs[k] = v; },
    getAttribute(k) { return this.attrs[k]; },
    querySelector() { return null; }};
buttons['pinsp-stage'] = stage;
const morphBar = {attrs: {'data-busy': '0'}, setAttribute(k, v) { this.attrs[k] = v; },
    getAttribute(k) { return this.attrs[k]; }};
buttons['pinsp-morph'] = morphBar;
// The partners button. Its `disabled` state is load-bearing -- a conformation
// with no partners of its own must not offer the toggle -- so the stub
// honours it the way a browser does and refuses the click.
const partnerBtn = stubButton();
partnerBtn.disabled = false;
partnerBtn.click = function () { if (!this.disabled) this.handlers.forEach((fn) => fn()); };
buttons['pinsp-partners'] = partnerBtn;

// One button per conformation, discovered the way the page discovers them.
const morphButtons = ((state.morph && state.morph.conformers) || []).map((c, i) => {
    const b = stubButton();
    b.attrs['data-conf'] = String(i);
    return b;
});

// The split grid: enough of a box model to exercise the clamp. The grid is
// 1000px wide starting at x=0, so a pointer at clientX=x asks for a panel of
// (1000 - x) px, which the page must clamp to [240, 1000-320].
const grid = {
    props: {'--pinsp-panel': '340px'},
    style: {setProperty(k, v) { grid.props[k] = v; }},
    getBoundingClientRect: () => ({left: 0, right: 1000, width: 1000, top: 0, bottom: 800, height: 800}),
};
const splitter = {
    attrs: {}, handlers: {},
    addEventListener(type, fn) { (this.handlers[type] = this.handlers[type] || []).push(fn); },
    setAttribute(k, v) { this.attrs[k] = v; },
    removeAttribute(k) { delete this.attrs[k]; },
    getAttribute(k) { return this.attrs[k]; },
    fire(type, ev) { (this.handlers[type] || []).forEach((fn) => fn(Object.assign(
        {preventDefault() {}, clientX: 0, shiftKey: false}, ev))); },
};
buttons['pinsp-split'] = splitter;
// ...and the sequence strip's own seam, which moves --pinsp-seq the way the
// vertical one moves --pinsp-panel. Same stub shape, a different axis: the
// grid's box is 800px tall and ends at y=800, so a pointer at clientY=y asks
// for a strip of (800 - y) px.
const hsplitter = {
    attrs: {}, handlers: {},
    addEventListener(type, fn) { (this.handlers[type] = this.handlers[type] || []).push(fn); },
    setAttribute(k, v) { this.attrs[k] = v; },
    removeAttribute(k) { delete this.attrs[k]; },
    getAttribute(k) { return this.attrs[k]; },
    fire(type, ev) { (this.handlers[type] || []).forEach((fn) => fn(Object.assign(
        {preventDefault() {}, clientY: 0, shiftKey: false}, ev))); },
};
buttons['pinsp-hsplit'] = hsplitter;
// BOTH ROWS THE SEAM SPLITS. The page seeds these in its own stylesheet from
// the caller's `display` height, so the stub reads the same number out of the
// state rather than inventing one -- otherwise the arithmetic under test and
// the arithmetic in the harness are two different sums.
grid.props['--pinsp-seq'] = '152px';
grid.props['--pinsp-view'] = String(
    (((state.viewer || {}).config || {}).display || {}).size &&
    state.viewer.config.display.size[1] || 720) + 'px';

// THE STRIP'S BODY, WHICH READS BACK THE MARKUP THE PAGE WROTE.
//
// The page builds the strip as one innerHTML string -- 20,000 createElement
// calls is not a thing to do on load -- so the honest stub is one that parses
// what it was handed rather than one that pretends to be a DOM. Cells come
// back from querySelectorAll in the order they were written, carrying their
// letter, their author-number tick and their title, so a test can assert what
// the reader will actually see.
function stubCell(index, attrs, text) {
    return {
        tag: 'span', style: {}, attrs: Object.assign({'data-i': String(index)}, attrs),
        textContent: text,
        setAttribute(k, v) { this.attrs[k] = v; },
        removeAttribute(k) { delete this.attrs[k]; },
        getAttribute(k) { return this.attrs[k] === undefined ? null : this.attrs[k]; },
    };
}
const seqBody = {
    html: '', cells: [], labels: [], handlers: {},
    set innerHTML(value) {
        this.html = value;
        this.cells = [];
        this.labels = [];
        const cell = /<span class="pinsp-cell" data-i="(\d+)"([^>]*)>([^<]*)<\/span>/g;
        let m;
        while ((m = cell.exec(value)) !== null) {
            const attrs = {};
            const extra = /([a-z-]+)="([^"]*)"/g;
            let a;
            while ((a = extra.exec(m[2])) !== null) attrs[a[1]] = a[2];
            this.cells.push(stubCell(Number(m[1]), attrs, m[3]));
        }
        const label = /<button type="button" class="pinsp-seqlab" data-run="(\d+)"[^>]*>([^<]*)</g;
        while ((m = label.exec(value)) !== null) {
            this.labels.push(stubCell(-1, {'data-run': m[1], 'data-i': undefined}, m[2]));
        }
        for (const l of this.labels) delete l.attrs['data-i'];
    },
    get innerHTML() { return this.html; },
    querySelectorAll(sel) { return sel === '.pinsp-cell' ? this.cells : []; },
    addEventListener(type, fn) { (this.handlers[type] = this.handlers[type] || []).push(fn); },
    fire(type, ev) { (this.handlers[type] || []).forEach((fn) => fn(Object.assign(
        {preventDefault() {}, shiftKey: false, metaKey: false, ctrlKey: false}, ev))); },
};
buttons['pinsp-seqbody'] = seqBody;
// The bar: a readout, a distance box and three buttons. `disabled` is
// load-bearing on all three -- a tool that acts on the selection must not be
// pressable when there is none -- so the stubs honour it as a browser does.
const selOut = {textContent: ''};
buttons['pinsp-sel'] = selOut;
// ONE `sel` BUTTON PER LAYER ROW, discovered by attribute the way the page
// discovers them. `shiftClick` is a separate entry point because the handler
// reads the modifier off the event, and that is the whole of how several
// layers are selected together.
const layerPicks = state.layers.map((layer) => {
    const b = stubButton();
    b.attrs['data-select'] = layer.layer_id;
    b.click = function () { this.handlers.forEach((fn) => fn(
        {preventDefault() {}, stopPropagation() {}, shiftKey: false})); };
    b.shiftClick = function () { this.handlers.forEach((fn) => fn(
        {preventDefault() {}, stopPropagation() {}, shiftKey: true})); };
    return b;
});
const aroundInput = {value: '5'};
buttons['pinsp-around-a'] = aroundInput;
for (const id of ['pinsp-around', 'pinsp-sc-on', 'pinsp-sc-off', 'pinsp-sel-clear',
                  'pinsp-layers-select']) {
    const b = stubButton();
    b.disabled = false;
    b.click = function () { if (!this.disabled) this.handlers.forEach((fn) => fn()); };
    buttons[id] = b;
}
global.getComputedStyle = (el) => ({getPropertyValue: (k) => (el.props || {})[k] || ''});

function bySelector(sel) {
    if (sel === '.protein-inspector') return grid;
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
// REAL REGISTRIES, NOT NO-OPS. Both of these were `addEventListener() {}`,
// which let every document-level handler the page installs pass untested --
// and the selection-change listener is now the single place the strip, the
// readout and the three tool buttons are refreshed from, so a harness that
// cannot deliver that event cannot see any of it. The window registry matters
// for the same reason: a drag in the strip ends on a window mouseup.
const winHandlers = {};
const docHandlers = {};
function fireWindow(type, ev) {
    (winHandlers[type] || []).forEach((fn) => fn(Object.assign(
        {preventDefault() {}, clientX: 0, clientY: 0, shiftKey: false}, ev)));
}
function fireDoc(type, ev) {
    (docHandlers[type] || []).forEach((fn) => fn(Object.assign(
        {preventDefault() {}, stopPropagation() {}}, ev)));
}
// BOTH PAINTERS, as the notebook bundle ships them -- it is the one build
// carrying the pair, which is what makes the painter a real choice there and
// not a request for a file that is absent. The page only offers the switch
// when both are present AND WebGL2 answers, so the stub has to answer both
// questions; `invalidate` is counted, because switching painters without
// dropping the mesh re-draws stale geometry.
let gpuInvalidated = 0;
global.window = {
    addEventListener(type, fn) { (winHandlers[type] = winHandlers[type] || []).push(fn); },
    removeEventListener() {},
    innerHeight: spec.innerHeight || 900,
    py2dmolCartoonGPU: {available: () => true, invalidate() { gpuInvalidated += 1; }},
    py2dmolCartoonPaint: function () {},
    performance: {now: () => clock}};
global.requestAnimationFrame = global.requestAnimationFrame || (fn => setTimeout(fn, 0));
global.CSS = {escape: (s) => s};
global.requestAnimationFrame = (fn) => rafQueue.push(fn);
// --- the viewer's own Style panel, as parts/panel.js builds it ------------
//
// Seven rows of `div.toggle-item`, each holding `div.half` cells, in the order
// STYLE_PANEL_ROWS declares them. The page folds this panel (foldStylePanel),
// and folding means moving row elements and one half-cell between parents --
// so the stub has to be a small real tree, not a bag of ids: a stub that
// answered querySelector without owning children would let a fold that moves
// nothing pass.
//
// Detail sits in the fifth row beside Pencil, Shadow, Ink and Ortho, inside
// its own `.half`, because that is the cell the page lifts out.
function stubNode(tag, className) {
    return {
        // `id` is a PROPERTY, not an attribute. The page assigns `node.id =`,
        // which is what a real DOM reflects into both; a stub that only kept
        // `attrs.id` matched the rows it had built itself and nothing the page
        // created, so a control the page added looked like a control missing.
        tag, className: className || '', attrs: {}, textContent: '', id: '',
        children: [], parentNode: null, title: '', type: '', checked: false,
        appendChild(node) {
            if (node.parentNode) {
                const at = node.parentNode.children.indexOf(node);
                if (at >= 0) node.parentNode.children.splice(at, 1);
            }
            node.parentNode = this;
            this.children.push(node);
            return node;
        },
        insertBefore(node, ref) {
            if (node.parentNode) {
                const at = node.parentNode.children.indexOf(node);
                if (at >= 0) node.parentNode.children.splice(at, 1);
            }
            node.parentNode = this;
            const where = this.children.indexOf(ref);
            this.children.splice(where < 0 ? this.children.length : where, 0, node);
            return node;
        },
        setAttribute(k, v) { this.attrs[k] = v; },
        getAttribute(k) { return this.attrs[k] === undefined ? null : this.attrs[k]; },
        handlers: {},
        addEventListener(type, fn) { (this.handlers[type] = this.handlers[type] || []).push(fn); },
        fire(type) { (this.handlers[type] || []).forEach((fn) => fn({target: this})); },
        // Depth-first over the subtree, matching only `#id` and `.class` --
        // the two forms foldStylePanel uses.
        querySelector(sel) {
            for (const kid of this.children) {
                if (sel[0] === '#' && kid.id === sel.slice(1)) return kid;
                if (sel[0] === '.' && String(kid.className).split(' ').includes(sel.slice(1))) {
                    return kid;
                }
                const deep = kid.querySelector(sel);
                if (deep) return deep;
            }
            return null;
        },
    };
}
const stylePanel = stubNode('div', '');
stylePanel.id = 'stylePanel';
const STYLE_ROWS = [
    ['styleSelect'],
    ['lineWidthSlider', 'outlineWidthSlider'],
    ['thicknessSlider', 'sheetFlatSlider'],
    ['highlightSlider', 'shadeSlider'],
    ['pencilSlider', 'shadowSlider', 'outlineTintSlider', 'detailSlider', 'orthoSlider'],
    ['colorSelect', 'selectionMarkSelect'],
    ['cartoonArrowsToggle', 'basePlatesToggle'],
];
for (const ids of STYLE_ROWS) {
    const row = stylePanel.appendChild(stubNode('div', 'toggle-item'));
    for (const id of ids) {
        const half = row.appendChild(stubNode('div', 'half'));
        half.appendChild(stubNode('input', '')).id = id;
    }
}
// WHAT THE READER IS LEFT LOOKING AT, read back off the tree: the ids still
// in a direct child of the panel, in order. Anything inside the fold is one
// click away and does not count.
const styleFront = () => stylePanel.children
    .filter((row) => row.className !== 'pinsp-fine')
    .map((row) => row.children.map((cell) => (cell.children[0] || {}).id))
    .reduce((all, ids) => all.concat(ids), [])
    .filter((id) => !!id);
buttons.stylePanel = stylePanel;

global.document = {
    createElement: (tag) => stubNode(tag, ''),
    getElementById: (id) => (id === 'pinsp-inspection-state'
        ? {textContent: stateMatch[1]}
        : (id === 'pinsp-residue-details' ? detailsNode
            // ...and the viewer's own controls, which live in the Style panel
            // tree rather than in `buttons`. getElementById is document-wide
            // in a browser, so a stub that only knew the page's own ids made
            // every lookup of a VIEWER control return null -- and a listener
            // the page attaches to one of those silently never attached.
            : (id === 'pinsp-card' ? cardNode
                : (buttons[id] || stylePanel.querySelector(`#${id}`) || null)))),
    querySelector: bySelector,
    querySelectorAll: (sel) => (sel === '[data-annotation]'
        ? Object.values(annotationRows)
        : (sel === '[data-layer]' ? Object.values(layerBoxes)
            : (sel === '.bm-btn[data-conf]' ? morphButtons
                : (sel === '[data-select]' ? layerPicks
                    : (sel === '.pinsp-cell' ? seqBody.cells : []))))),
    addEventListener(type, fn) { (docHandlers[type] = docHandlers[type] || []).push(fn); },
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
    // PICKING IS OFF IN THE NOTEBOOK BUNDLE (src/parts/ui.js turns it on for
    // Focus mode only), so an exported page has to switch it on itself. False
    // here, as it is there.
    selectionEnabled: false,
    positionNames: spec.positionNames || chains.map(() => 'ALA'),
    positionTypes: spec.positionTypes || chains.map(() => 'P'),
    localRangeOf() { return {off: 0, end: chains.length}; },
    // py2Dmol's own setter dispatches on DOCUMENT and that is how every
    // surface -- the canvas, the strip, the tools, the manifest rows -- stays
    // in step. A stub that only assigned the field tested none of it.
    setResidueSelection(next) {
        const many = next && (next.size || (Array.isArray(next) && next.length));
        this.residueSelection = many ? new Set(next) : new Set();
        fireDoc('py2dmol-residue-selection-change');
    },
    clearResidueSelection() {
        this.residueSelection = new Set();
        fireDoc('py2dmol-residue-selection-change');
    },
    // A NEIGHBOURHOOD, FAKED BY SEQUENCE DISTANCE. The real one measures atom
    // to atom on a spatial grid; what the page has to get right is that it
    // passes the seed and the cutoff and adopts whatever comes back, so the
    // stub answers something checkable instead: one position per angstrom
    // either side, seed included, as residuesWithin includes it.
    withinCalls: [],
    residuesWithin(seed, cutoff) {
        this.withinCalls.push({seed: Array.from(seed).sort((a, b) => a - b), cutoff});
        const out = new Set(seed);
        const reach = Math.round(Number(cutoff));
        for (const i of seed) {
            for (let d = -reach; d <= reach; d++) {
                const k = i + d;
                if (k >= 0 && k < chains.length) out.add(k);
            }
        }
        return out;
    },
    sidechainCalls: [],
    // The renderer's own merged answer, which the page asks for first.
    shownSidechainSet() {
        return (object.sidechains instanceof Set) ? object.sidechains : new Set();
    },
    showSidechains(sel) {
        this.sidechainCalls.push({on: true, positions: (sel.positions || []).slice()});
        const cur = (object.sidechains instanceof Set) ? object.sidechains : new Set();
        for (const i of sel.positions || []) cur.add(i);
        object.sidechains = cur;
        return this;
    },
    hideSidechains(sel) {
        this.sidechainCalls.push({on: false, positions: (sel.positions || []).slice()});
        const cur = (object.sidechains instanceof Set) ? object.sidechains : new Set();
        for (const i of sel.positions || []) cur.delete(i);
        object.sidechains = cur;
        return this;
    },
    visibility: null,
    showAll() { this.visibility = 'all'; },
    chainKeyAt(i) { return 'obj:' + this.chains[i]; },
    // py2Dmol's own setVisibility sends the viewer back to the first frame.
    // The stub reproduces that so the page's restore is actually exercised;
    // without it the regression is invisible here and visible to the reader.
    // It drops the frame inline AND again on the next animation frame: a
    // restore that only runs synchronously loses to the late one, which is
    // what the first fix did.
    setVisibility(patch) {
        this.visibility = patch;
        this.currentFrame = 0;
        requestAnimationFrame(() => { this.currentFrame = 0; });
    },
    setFrame(i) { this.currentFrame = i; this.framesVisited.push(i);
        if (this.onFrame) this.onFrame(i); renders++; },
    updateUIControls() { uiUpdates++; },
    // WHICH PAINTER, as the bundle's config seeds it (display gpu, default
    // True). The page's painter switch assigns this, so a stub without it
    // would let a switch that sets nothing pass.
    useGPU: true,
    // CARTOON SAMPLING, AND WHAT A STYLE SWITCH DOES TO IT. core/mol.js's
    // _applyLookDefaults re-asserts every style-owned control from
    // LOOK_DEFAULTS, and all three entries there carry detail 4 -- with no
    // "did a person choose this" latch of the kind width and thickness have.
    // So setStyle stamping 4 here is not a stub simplification, it is the
    // behaviour the page has to survive: a bundle opens as tube, and the
    // reader's first move is to pick Richardson.
    cartoonDetail: (((state.viewer || {}).config || {}).rendering || {}).detail,
    styleChanges: [],
    setStyle(name) {
        this.style = name === 'tube' ? 'tube' : 'cartoon';
        this.stylePreset = name === 'tube' ? null : name;
        this.cartoonDetail = 4;
        this.styleChanges.push(name);
        renders++;
    },
    setPreset(name) { this.setStyle(name); },
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
    report.morph = {
        framesInFile: state.morph.conformers.length,
        uiUpdates,
        pressedBefore: morphButtons.map((b) => b.getAttribute('aria-pressed')),
    };
    renderer.framesVisited.length = 0;
    // Record what the buffer held at each drawn step, so "went straight there"
    // can be told apart from "walked through the middle conformation".
    const trail = [];
    renderer.onFrame = (i) => {
        const f = object.frames[i];
        if (f) trail.push({index: i, y: f.coords[0][1]});
    };
    morphButtons[target].click();
    drain(400);
    renderer.onFrame = null;
    const visited = renderer.framesVisited;
    report.morph.framesAfterExpansion = object.frames.length;
    report.morph.visitedCount = visited.length;
    report.morph.landedOn = visited[visited.length - 1];
    report.morph.expectedLanding = target;
    // Every frame drawn DURING the animation is the scratch buffer; only the
    // final frame is a conformation. Passing through frame 1 would mean the
    // old chain behaviour is back.
    report.morph.buffersOnly = visited.slice(0, -1)
        .every((v) => v === state.morph.conformers.length);
    report.morph.bufferMaxY = Math.max(...trail.map((t) => t.y));
    report.morph.bufferEndY = trail[trail.length - 1].y;
    report.morph.pressedAfter = morphButtons.map((b) => b.getAttribute('aria-pressed'));
    report.morph.note = buttons['pinsp-morph-note'].textContent;
}

// --- partners --------------------------------------------------------------
// Pressing the button must hide exactly the current conformation's partner
// positions; and a partner must follow the conformation, not linger.
if (state.partners && state.partners.owner) {
    const ownerOf = state.partners.owner;
    // The morph section above left the viewer on the last conformation; come
    // back to the first so the partner assertions start from a known state.
    if (state.morph && morphButtons.length > 1) {
        morphButtons[0].click();
        drain(400);
    }
    const visible = () => {
        const v = renderer.visibility;
        return (v && v.positions) ? Array.from(v.positions).sort((a, b) => a - b) : null;
    };
    const targetOnly = chains.map((c, i) => [c, i])
        .filter(([c]) => ownerOf[c] === undefined).map(([, i]) => i);
    const ownedBy = (k) => chains.map((c, i) => [c, i])
        .filter(([c]) => ownerOf[c] === k).map(([, i]) => i);

    report.partners = {
        disabledAtStart: partnerBtn.disabled,
        pressedAtStart: partnerBtn.getAttribute('aria-pressed'),
        onAtStart: visible(),
        expectedOn: targetOnly.concat(ownedBy(0)).sort((a, b) => a - b),
    };
    renderer.framesVisited.length = 0;
    partnerBtn.click();
    report.partners.framesTouchedByToggle = renderer.framesVisited.slice();
    report.partners.frameAfterToggle = renderer.currentFrame;
    report.partners.afterHide = visible();
    report.partners.expectedOff = targetOnly;
    report.partners.pressedAfterHide = partnerBtn.getAttribute('aria-pressed');
    partnerBtn.click();
    report.partners.afterShow = visible();

    // ...and on a state that is not the first: toggling there must not send
    // the viewer home. This is the regression the reader actually hit. The
    // expected frame is the conformation's OWN index -- reading currentFrame
    // back after the morph would capture whatever the bug left behind and
    // compare it against itself.
    if (state.morph && morphButtons.length > 1) {
        const last = morphButtons.length - 1;
        morphButtons[last].click();
        drain(400);
        partnerBtn.disabled = false;
        partnerBtn.click();
        drain(6);
        report.partners.frameAfterHideOnState = renderer.currentFrame;
        partnerBtn.click();
        drain(6);
        report.partners.frameAfterShowOnState = renderer.currentFrame;
        report.partners.expectedFrameOnState = last;
        morphButtons[0].click();
        drain(400);
    }

    // Morph to the last conformation: its own partners appear, state 0's go.
    if (state.morph && morphButtons.length > 1) {
        const last = morphButtons.length - 1;
        morphButtons[last].click();
        drain(400);
        report.partners.afterMorph = visible();
        report.partners.expectedAfterMorph = targetOnly.concat(ownedBy(last))
            .sort((a, b) => a - b);
        report.partners.disabledAfterMorph = partnerBtn.disabled;
    }
}

// --- the split boundary ----------------------------------------------------
// Dragging the seam must move the panel and clamp at both ends, so the two
// panes can neither overlap nor leave the stage unusably narrow.
if (splitter.handlers.pointerdown) {
    const panel = () => parseFloat(grid.props['--pinsp-panel']);
    splitter.fire('pointerdown', {pointerId: 1});
    report.split = {dragAttr: splitter.getAttribute('data-drag')};
    splitter.fire('pointermove', {clientX: 600});
    report.split.at600 = panel();
    splitter.fire('pointermove', {clientX: 100});   // asks for 900, too wide
    report.split.clampedWide = panel();
    splitter.fire('pointermove', {clientX: 990});   // asks for 10, too narrow
    report.split.clampedNarrow = panel();
    splitter.fire('pointerup', {});
    // JSON.stringify drops undefined keys, so say null out loud.
    report.split.dragAttrAfter = splitter.getAttribute('data-drag') || null;
    // ...and a pointermove after release must not move it.
    splitter.fire('pointermove', {clientX: 500});
    report.split.afterRelease = panel();
    splitter.fire('keydown', {key: 'ArrowLeft'});
    report.split.afterArrowLeft = panel();
    splitter.fire('keydown', {key: 'Home'});
    report.split.afterHome = panel();
}

// --- Detail survives a style switch ---------------------------------------
// The whole point of rendering a bundle at detail 8: a reader who picks
// Richardson -- the style they would pick it FOR -- must not silently drop to
// the library default of 4.
{
    const want = (((state.viewer || {}).config || {}).rendering || {}).detail;
    report.detailKeeper = {configured: want, atStart: renderer.cartoonDetail};
    renderer.setStyle('richardson');
    drain(4);
    report.detailKeeper.afterRichardson = renderer.cartoonDetail;
    // ...and the reader still outranks the file: once the slider is dragged,
    // their number is the one that survives the next switch.
    const slider = stylePanel.querySelector('#detailSlider');
    if (slider) {
        slider.value = '3';
        slider.fire('input');
        renderer.setStyle('3d');
        drain(4);
        report.detailKeeper.afterDragThenSwitch = renderer.cartoonDetail;
    }
    renderer.setStyle('tube');
    drain(4);
}

// --- the folded Style panel -----------------------------------------------
// Four controls out front -- Style, Detail, Color, Sele -- and every other row
// inside one closed <details>. Asserted on the TREE rather than on a count,
// because the failure that matters is a control that went missing rather than
// one that moved: `front` plus `folded` must still be the whole panel.
{
    const fold = stylePanel.querySelector('.pinsp-fine');
    report.styleFold = {
        front: styleFront(),
        foldedRows: fold ? fold.children.filter((n) => n.tag !== 'summary').length : 0,
        folded: fold
            ? fold.children
                .filter((n) => n.tag !== 'summary')
                .map((row) => row.children.map((c) => (c.children[0] || {}).id))
                .reduce((all, ids) => all.concat(ids), [])
            : [],
        summary: fold ? (fold.children[0] || {}).textContent : null,
        closed: fold ? !fold.attrs.open : null,
        // ONE FOLD, however many frames the init loop ran. It calls
        // foldStylePanel on every animation frame until the viewer reports
        // ready, so a pass that did not notice its own previous work would
        // nest a fold inside a fold and empty the panel.
        foldCount: stylePanel.children.filter((n) => n.className === 'pinsp-fine').length,
    };
    // THE PAINTER SWITCH. Off must reach the renderer AND drop the GPU mesh:
    // switching without invalidating re-draws geometry built from state that
    // moved on while the painter was off.
    const gpuBox = stylePanel.querySelector('#pinsp-gpu');
    if (gpuBox) {
        report.styleFold.gpuBefore = renderer.useGPU;
        gpuBox.checked = false;
        gpuBox.fire('change');
        report.styleFold.gpuAfterOff = renderer.useGPU;
        report.styleFold.invalidated = gpuInvalidated;
        gpuBox.checked = true;
        gpuBox.fire('change');
        report.styleFold.gpuAfterOn = renderer.useGPU;
    }
    report.styleFold.detailIsFront = report.styleFold.front.includes('detailSlider');
}

// --- controls collapse ----------------------------------------------------
buttons['pinsp-controls'].click();
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
    buttons['pinsp-layers-none'].click();
    report.untickAll = {
        paint: paint(),
        checked: Object.values(layerBoxes).filter((b) => b.checked).length,
        hiddenRows: Object.entries(annotationRows).filter(([, r]) => r.hidden).length,
    };
    buttons['pinsp-layers-all'].click();
    report.tickAll = {
        paint: paint(),
        checked: Object.values(layerBoxes).filter((b) => b.checked).length,
    };
    const someAnnotation = state.annotations.find((a) => a.resolved);
    if (someAnnotation) window.proteinInspector.selectAnnotation(someAnnotation.annotation_id);
    report.afterSelect = {details: detailsNode.textContent.length,
                          selection: renderer.residueSelection.size};
    buttons['pinsp-clear-residue'].click();
    report.afterClear = {details: detailsNode.textContent,
                         selection: renderer.residueSelection.size};
}

// --- the sequence strip ---------------------------------------------------
// What the page WROTE: one cell per drawn position, in file order, carrying
// the letter, the author-number tick every tenth residue, and the layer
// colour the structure has.
const sel = () => Array.from(renderer.residueSelection).sort((a, b) => a - b);
const marked = () => seqBody.cells.filter((c) => c.getAttribute('data-sel'))
    .map((c) => Number(c.getAttribute('data-i'))).sort((a, b) => a - b);
report.seq = {
    cells: seqBody.cells.length,
    letters: seqBody.cells.map((c) => c.textContent).join(''),
    ticks: seqBody.cells.filter((c) => c.getAttribute('data-num'))
        .map((c) => c.getAttribute('data-num')),
    chainLabels: seqBody.labels.map((l) => l.textContent),
    firstTitle: seqBody.cells.length ? seqBody.cells[0].getAttribute('title') : null,
    colours: seqBody.cells.map((c) => c.style.background || null),
    // The flag core/mol.js gates every click on. False in the bundle as it
    // ships; the page has to turn it on or the canvas is inert.
    selectionEnabled: renderer.selectionEnabled,
};

if (seqBody.cells.length > 9) {
    // ...a click on a letter selects that residue, marks it, names it, and
    // opens the card.
    seqBody.fire('mousedown', {target: seqBody.cells[2]});
    report.seq.afterClick = {selection: sel(), marked: marked(),
        readout: selOut.textContent, card: cardNode.innerHTML, cardHidden: cardNode.hidden};
    // ...dragging to another letter takes the range between them.
    seqBody.fire('mousemove', {target: seqBody.cells[5]});
    report.seq.afterDrag = sel();
    // ...and a release ends the drag, so moving over a sixth letter is hover
    // and not selection.
    fireWindow('mouseup', {});
    seqBody.fire('mousemove', {target: seqBody.cells[8]});
    report.seq.afterRelease = sel();
    // ...A CLICK NEVER REPLACES. Pressing an unselected letter ADDS it to
    // whatever was there; pressing a selected one takes it away. One rule,
    // both directions, no modifier.
    seqBody.fire('mousedown', {target: seqBody.cells[8]});
    report.seq.afterAddClick = sel();
    seqBody.fire('mousedown', {target: seqBody.cells[8]});
    report.seq.afterRemoveClick = sel();
    // ...and a drag STARTED ON A SELECTED letter subtracts its whole range,
    // leaving everything outside the range alone.
    seqBody.fire('mousedown', {target: seqBody.cells[3]});
    report.seq.beforeSubtract = sel();
    seqBody.fire('mousemove', {target: seqBody.cells[5]});
    report.seq.afterSubtractDrag = sel();
    // ...dragging BACK over your own path shrinks the range rather than
    // ratcheting it: each move recomputes from the mousedown snapshot.
    seqBody.fire('mousemove', {target: seqBody.cells[4]});
    report.seq.afterDragBack = sel();
    fireWindow('mouseup', {});
    // ...shift extends from the last letter pressed.
    seqBody.fire('mousedown', {target: seqBody.cells[1]});
    seqBody.fire('mousedown', {target: seqBody.cells[4], shiftKey: true});
    report.seq.afterShift = sel();
    // ...and the chain label takes the whole chain.
    if (seqBody.labels.length) {
        seqBody.fire('mousedown', {target: seqBody.labels[0]});
        report.seq.afterChainLabel = renderer.residueSelection.size;
        // ...and pressing it again takes the whole chain back out.
        seqBody.fire('mousedown', {target: seqBody.labels[0]});
        report.seq.afterChainLabelAgain = renderer.residueSelection.size;
    }

    // --- select around ----------------------------------------------------
    // The distance box is read, the renderer's own search is asked, and
    // whatever it answers becomes the selection.
    seqBody.fire('mousedown', {target: seqBody.cells[4]});
    renderer.withinCalls.length = 0;
    aroundInput.value = '2';
    report.around = {disabledWithSelection: buttons['pinsp-around'].disabled};
    buttons['pinsp-around'].click();
    report.around.calls = renderer.withinCalls.slice();
    report.around.selection = sel();
    report.around.readout = selOut.textContent;
    // ...and with nothing selected it is not pressable at all.
    buttons['pinsp-sel-clear'].click();
    report.around.disabledWhenEmpty = buttons['pinsp-around'].disabled;
    report.around.clearedSelection = sel();

    // --- side chains of the selection -------------------------------------
    buttons['pinsp-sel-clear'].click();
    seqBody.fire('mousedown', {target: seqBody.cells[6]});
    renderer.sidechainCalls.length = 0;
    const scOn = buttons['pinsp-sc-on'], scOff = buttons['pinsp-sc-off'];
    const pressed = () => [scOn.getAttribute('aria-pressed'),
                           scOff.getAttribute('aria-pressed')];
    report.sidechains = {
        stateSaysPresent: !!state.sidechains,
        disabled: scOn.disabled,
        pressed: pressed(),
    };
    scOn.click();
    report.sidechains.afterShow = {calls: renderer.sidechainCalls.slice(),
        pressed: pressed()};
    scOff.click();
    report.sidechains.afterHide = {calls: renderer.sidechainCalls.slice(),
        pressed: pressed()};
    // ...THE THIRD STATE. Draw one residue's side chain, then select it
    // together with one that has none: the pair must fill NEITHER button,
    // which is the state a single latch had to report as "off".
    scOn.click();
    seqBody.fire('mousedown', {target: seqBody.cells[7]});
    report.sidechains.mixedPressed = pressed();
    report.sidechains.mixedSelection = sel();
    // ...and Show resolves a disagreeing selection in ONE press rather than
    // flipping it to whichever state the majority was not in.
    renderer.sidechainCalls.length = 0;
    scOn.click();
    report.sidechains.afterResolve = {calls: renderer.sidechainCalls.slice(),
        pressed: pressed()};
    // ...and the per-residue state outlives the selection: dropping it must
    // not undraw anything, and the readout has to say so.
    buttons['pinsp-sel-clear'].click();
    report.sidechains.afterDeselect = {
        shown: renderer.shownSidechainSet ? renderer.shownSidechainSet().size : null,
        readout: selOut.textContent, disabled: scOn.disabled, pressed: pressed()};

    // --- a click in the VIEWER -------------------------------------------
    // core/mol.js's mouseup handler calls setResidueSelection and nothing
    // else; everything the reader sees has to follow from the event.
    // ...INCLUDING THE REPAINT. mol.js's mouseup renders only for a large
    // molecule, so on anything smaller the mark appeared on the first frame
    // after the click -- nudge the structure and there it was.
    const rendersBefore = renders;
    renderer.setResidueSelection(new Set([7]));
    report.viewerClick = {card: cardNode.innerHTML, details: detailsNode.textContent,
        marked: marked(), readout: selOut.textContent,
        rendered: renders - rendersBefore};
    // ...a whole chain (double-click) reports the set rather than a card.
    renderer.setResidueSelection(new Set([0, 1, 2]));
    report.viewerClick.chainCard = cardNode.hidden;
    report.viewerClick.chainReadout = selOut.textContent;
    // ...and a click on the background clears all of it.
    renderer.clearResidueSelection();
    report.viewerClick.afterBackground = {cardHidden: cardNode.hidden,
        details: detailsNode.textContent, marked: marked(),
        readout: selOut.textContent};

    // --- a layer AS a selection -------------------------------------------
    // The tools act on a selection, and a layer is already a named set of
    // residues -- so the two halves of the page were not connected. `sel`
    // selects one layer, shift adds, and the header's Select takes every
    // ticked layer at once, which is what "side chains for these" means.
    buttons['pinsp-sel-clear'].click();
    if (layerPicks.length) {
        layerPicks[0].click();
        report.layerSelect = {one: sel(), readout: selOut.textContent,
            marked: marked()};
        layerPicks[layerPicks.length - 1].click();
        report.layerSelect.replaced = sel();
        layerPicks[0].shiftClick();
        report.layerSelect.added = sel();
    }
    if (buttons['pinsp-layers-select'] && state.layers.length) {
        buttons['pinsp-sel-clear'].click();
        buttons['pinsp-layers-select'].click();
        report.layerSelect = report.layerSelect || {};
        report.layerSelect.allVisible = sel();
        // ...and with a layer unticked, its residues are not in the answer.
        const firstLayer = state.layers[0].layer_id;
        layerBoxes[firstLayer].checked = false;
        layerBoxes[firstLayer].handlers.forEach((fn) => fn());
        buttons['pinsp-layers-select'].click();
        report.layerSelect.visibleAfterUntick = sel();
        layerBoxes[firstLayer].checked = true;
        layerBoxes[firstLayer].handlers.forEach((fn) => fn());
        // ...and the side-chain button then acts on the whole set.
        buttons['pinsp-layers-select'].click();
        renderer.sidechainCalls.length = 0;
        buttons['pinsp-sc-on'].click();
        report.layerSelect.sidechainsForLayers =
            (renderer.sidechainCalls[0] || {}).positions || null;
        buttons['pinsp-sc-off'].click();
    }

    // --- the hover bridge -------------------------------------------------
    // core/mol.js calls window.SEQ.setHoveredResidue while the pointer is
    // over the structure, when a strip is present to receive it.
    if (window.SEQ && window.SEQ.setHoveredResidue) {
        const num = residueNumbers[5];
        window.SEQ.setHoveredResidue({chain: chains[5], resSeq: num, resName: 'ALA'});
        report.hover = {lit: seqBody.cells.filter((c) => c.getAttribute('data-hover'))
            .map((c) => Number(c.getAttribute('data-i')))};
        window.SEQ.setHoveredResidue(null);
        report.hover.afterLeave = seqBody.cells
            .filter((c) => c.getAttribute('data-hover')).length;
    }
}

// --- the sequence strip's own boundary -------------------------------------
// ONE SEAM, TWO ROWS, CONSTANT SUM. What is reported at every step is the PAIR,
// because the bug was that only one of them moved: the strip grew, the page
// got taller, and the only thing the reader saw move was the credit line.
//
// The gesture is anchored to pointerdown. The first version measured against
// the grid's own bottom edge, which moves when the strip resizes, so each
// event read back its own previous result and the seam ran away from the
// pointer. A harness whose grid rect is a constant cannot see that -- so the
// travel asserted below is deliberately pointer-relative, and `atAnchor`
// pins the no-movement case that a feedback loop could never hold still.
if (hsplitter.handlers.pointerdown) {
    const pair = () => ({seq: parseFloat(grid.props['--pinsp-seq']),
        view: parseFloat(grid.props['--pinsp-view'])});
    report.hsplit = {start: pair()};
    hsplitter.fire('pointerdown', {pointerId: 2, clientY: 500});
    report.hsplit.dragAttr = hsplitter.getAttribute('data-drag');
    // Still at the anchor: nothing has moved, and nothing drifts.
    hsplitter.fire('pointermove', {clientY: 500});
    report.hsplit.atAnchor = pair();
    // 40px up: the strip takes 40, the picture gives up 40.
    hsplitter.fire('pointermove', {clientY: 460});
    report.hsplit.up40 = pair();
    // ...and back to the anchor, from the anchor, not from where it got to.
    hsplitter.fire('pointermove', {clientY: 500});
    report.hsplit.backToAnchor = pair();
    // Far down: the strip clamps at 64 and the picture keeps the rest.
    hsplitter.fire('pointermove', {clientY: 1200});
    report.hsplit.clampedShort = pair();
    // Far up: the picture clamps at 200 and the strip keeps the rest.
    hsplitter.fire('pointermove', {clientY: -400});
    report.hsplit.clampedTall = pair();
    hsplitter.fire('pointerup', {});
    report.hsplit.dragAttrAfter = hsplitter.getAttribute('data-drag') || null;
    hsplitter.fire('pointermove', {clientY: 300});
    report.hsplit.afterRelease = pair();
    // Home first, then a step: pressing a step key while the strip is already
    // clamped at its tallest asserts nothing about the step.
    hsplitter.fire('keydown', {key: 'Home'});
    report.hsplit.afterHome = pair();
    hsplitter.fire('keydown', {key: 'ArrowUp'});
    report.hsplit.afterArrowUp = pair();
    hsplitter.fire('keydown', {key: 'ArrowDown'});
    report.hsplit.afterArrowDown = pair();
}

report.finalRenders = renders;
process.stdout.write(JSON.stringify(report));
