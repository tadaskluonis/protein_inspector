/* py2Dmol plugin: hydrogen bonds as dashed sticks.
 *
 * The payload is a list of bonds with explicit donor/acceptor coordinates in
 * model space. Whether they are drawn is ONE question with one answer for the
 * whole page - the host's toggle, read from window.__pinspHbonds - rather than
 * something the selection or the focus decides for the reader.
 */
(function () {
    var P = window.py2dmolPlugins = window.py2dmolPlugins || { list: [], pending: [] };

    // FOCUS OUTRANKS THE BUTTON, and only while it is on. The host publishes the
    // set it is focused on; these bonds are the reason the reader moved in, so
    // they are drawn whatever the button says, and only for that set. With no
    // focus the button is the whole answer.
    function focusSet() {
        var f = window.__pinspHbondsFocus;
        return (f && f.active && f.set && f.set.length) ? f.set : null;
    }

    function shownNow(ctx) {
        if (focusSet()) return true;
        if (typeof window.__pinspHbonds === 'boolean') return window.__pinspHbonds;
        var o = ctx && ctx.options;
        return !(o && o.shown === false);
    }

    // (chain, author number) -> position index, from the renderer's own arrays,
    // so a payload written in author numbering can be matched to a focused set
    // written in position indices.
    // 🔴 THE RESIDUE, NOT ONE OF ITS ATOMS. Showing side chains APPENDS their
    // atoms as real positions, and each appended atom carries its OWNER's chain
    // and author number in these same arrays - so a map written straight across
    // them ends on the last atom of the residue instead of the residue: A7
    // resolved to 106 on a 76-residue ubiquitin, the focused set is residue
    // indices, nothing matched, and the dashes silently did not draw. Skip the
    // appended atoms (sidechainMap knows them) and let the first writer win.
    function indexOf(ctx) {
        var r = ctx.renderer, map = {}, ch = r.chains || [], num = r.residueNumbers || [];
        var sm = (r.sidechainMap && r.sidechainMap.has) ? r.sidechainMap : null;
        for (var i = 0; i < num.length; i++) {
            if (sm && sm.has(i)) continue;
            var key = String(ch[i]) + ':' + num[i];
            if (!(key in map)) map[key] = i;
        }
        return map;
    }

    function bondsOf(ctx) {
        var out = [];
        (ctx.payloads || []).forEach(function (p) {
            var b = p && p.payload && p.payload.bonds;
            if (b && b.length) out = out.concat(b);
        });
        return out;
    }

    function dist(a, b) {
        var dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    function lerp(a, b, t) {
        return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
    }

    var def = {
        name: 'hbonds',
        version: '3',
        apiVersion: 1,
        options: { shown: true, dash: 0.45, gap: 0.35, width: 0.13, color: '#2f7fb8', colors: null },

        init: function () { return {}; },

        // The toggle and the payload are all that change what is emitted.
        key: function (ctx) {
            var f = focusSet();
            return (shownNow(ctx) ? 'on' : 'off') + '|' + bondsOf(ctx).length + '|' + ctx.frame
                + '|' + (f ? f.slice().sort(function (a, b) { return a - b; }).join(',') : '');
        },

        prims: function (ctx) {
            var o = ctx.options || {}, bonds = bondsOf(ctx);
            // A COUNT THE PAGE CAN READ BACK, so a test (or a reader reporting a
            // fault) can say how many dashes the last frame actually drew.
            var focus = focusSet();
            window.__pinspHbondsLast = { shown: shownNow(ctx), focused: !!focus, bonds: 0, prims: 0 };
            if (!bonds.length || !shownNow(ctx)) return;
            if (focus) {
                var map = indexOf(ctx), want = {};
                for (var f = 0; f < focus.length; f++) want[focus[f]] = 1;
                bonds = bonds.filter(function (b) {
                    return want[map[b.ac + ':' + b.ar]] || want[map[b.bc + ':' + b.br]];
                });
                // WHY A FRAME DREW NOTHING, in the frame's own report: a focused
                // set that maps onto no residue is a mapping fault, and without
                // these three numbers it is indistinguishable from "this
                // neighbourhood genuinely has no hydrogen bonds".
                window.__pinspHbondsLast.mapped = Object.keys(map).length;
                window.__pinspHbondsLast.wanted = focus.length;
                window.__pinspHbondsLast.matched = bonds.length;
                window.__pinspHbondsLast.sample = (ctx.payloads[0] && ctx.payloads[0].payload
                    && ctx.payloads[0].payload.bonds[0])
                    ? [ctx.payloads[0].payload.bonds[0].ac + ':' + ctx.payloads[0].payload.bonds[0].ar,
                       map[ctx.payloads[0].payload.bonds[0].ac + ':' + ctx.payloads[0].payload.bonds[0].ar]]
                    : null;
                if (!bonds.length) return;
            }

            var dash = Math.max(0.1, o.dash || 0.45), gap = Math.max(0.05, o.gap || 0.35);
            var period = dash + gap, budget = (ctx.maxPrims || 8000) - 8, est = 0;
            for (var k = 0; k < bonds.length; k++) {
                est += Math.max(1, Math.round(dist(bonds[k].a, bonds[k].b) / period));
            }
            // Degrade to solid lines rather than let the budget drop the lot.
            var solid = est > budget, emitted = 0;
            var colors = o.colors || (ctx.payloads[0] && ctx.payloads[0].payload
                && ctx.payloads[0].payload.colors) || null;

            for (var j = 0; j < bonds.length; j++) {
                var bond = bonds[j];
                var col = (colors && colors[bond.k]) || o.color || '#2f7fb8';
                if (solid) { ctx.line(bond.a, bond.b, { color: col, width: o.width }); emitted++; continue; }
                var L = dist(bond.a, bond.b), n = Math.max(1, Math.round(L / period));
                var used = n * dash + (n - 1) * gap, start = Math.max(0, (L - used) / 2);
                for (var m = 0; m < n; m++) {
                    var t0 = (start + m * period) / L, t1 = (start + m * period + dash) / L;
                    ctx.line(lerp(bond.a, bond.b, t0), lerp(bond.a, bond.b, Math.min(1, t1)),
                        { color: col, width: o.width });
                    emitted++;
                }
            }
            window.__pinspHbondsLast = { shown: true, focused: !!focus, bonds: bonds.length, prims: emitted };
        },
    };

    if (P.register) P.register(def); else P.pending.push(def);
})();
