/* SOG1 comparative genomics report. All numbers come from results.js, written by the pipeline. */
(function () {
  "use strict";
  const D = window.SOG1;
  const $ = (s, r = document) => r.querySelector(s);
  const NS = "http://www.w3.org/2000/svg";
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const f = (x, d = 2) => (x == null || !isFinite(x) ? "–" : Number(x).toFixed(d));
  const pfmt = (p) => (p == null ? "–" : p < 0.001 ? "<0.001" : p.toFixed(3));
  const arr = (x) => (Array.isArray(x) ? x : x == null ? [] : [x]);

  if (!D) {
    const b = $("#banner"); b.hidden = false; b.textContent = "No results yet: the pipeline has not run. See the README.";
    return;
  }

  const LIN = { "Green alga": "--l-alga", "Charophyte alga": "--l-alga", Liverwort: "--l-bryo", Moss: "--l-bryo", Lycophyte: "--l-lyco",
    "Early angiosperm": "--l-early", Monocot: "--l-mono", Eudicot: "--l-dicot" };
  const linColor = (l) => css(LIN[l] || "--muted");
  const species = D.meta.species.filter((s) => s.included);
  const spBy = Object.fromEntries(D.meta.species.map((s) => [s.species, s]));
  const shortBy = Object.fromEntries(D.meta.species.map((s) => [s.short, s]));
  const clade = D.clade;
  const byKey = Object.fromEntries(clade.map((c) => [c.key, c]));
  const focusShort = spBy[D.meta.focus_species].short;
  const focus = clade.filter((c) => c.short === focusShort);
  const counts = Object.fromEntries(D.counts.map((c) => [c.species, c]));
  const LABELS = ["SOG1", "NAC044", "NAC085"];
  const labColor = { SOG1: css("--c1"), NAC044: css("--c2"), NAC085: css("--c3") };

  function svg(w, h) {
    const s = document.createElementNS(NS, "svg");
    s.setAttribute("viewBox", `0 0 ${w} ${h}`); s.setAttribute("width", w); s.setAttribute("height", h);
    s.setAttribute("role", "img");
    return s;
  }
  function add(parent, tag, attrs = {}, text) {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text != null) e.textContent = text;
    parent.appendChild(e);
    return e;
  }
  function table(el, head, rows) {
    el.innerHTML = "<thead><tr>" + head.map((h) => `<th class="${h.num ? "num" : ""}">${h.t}</th>`).join("") + "</tr></thead><tbody>" +
      rows.map((r) => "<tr>" + r.map((c, i) => `<td class="${head[i].num ? "num" : ""}">${c}</td>`).join("") + "</tr>").join("") + "</tbody>";
  }
  const geneName = (c) => `${c.short} ${c.gene}`;
  const set = (k, v) => document.querySelectorAll(`[data-k="${k}"]`).forEach((e) => (e.textContent = v));

  $("#generated").textContent = `Ensembl Plants release ${D.meta.release || "current"} · data downloaded ${D.meta.downloaded}`;
  set("release", D.meta.release || "current");
  set("nspecies", species.length);

  // ------------------------------------------------------------ headline numbers
  const fc = counts[D.meta.focus_species] || {};
  const sqRows = clade.filter((c) => c.sq_tail != null);
  const sqSig = sqRows.filter((c) => c.sq_tail_q < 0.05).length;
  const scd = sqRows.filter((c) => c.scd).length;
  const P = D.promoters.species || {};
  const M0 = D.promoters.motifs[0];
  const pMa = P[D.meta.focus_species] && P[D.meta.focus_species][M0];
  const sel = D.selection;
  const focusLab = LABELS.map((l) => [l, focus.filter((c) => c.clade_label === l).length]).filter((x) => x[1]);

  const kpis = [
    [species.length, "genomes searched", `${D.meta.species.length - species.length} not in this Ensembl release`],
    [fc.nac, "NAC genes in banana", fc.published ? `published: ${fc.published}` : ""],
    [focus.length, "SOG1-clade genes in banana", focusLab.map(([l, n]) => `${n} ${l}-like`).join(", ")],
    [`${sqSig}/${sqRows.length}`, "tails with more SQ/TQ than chance", `${scd} carry an SQ/TQ cluster`],
    [sel && sel.median_omega_domain != null ? `${f(sel.median_omega_domain)} vs ${f(sel.median_omega_tail)}` : "–", "median dN/dS, domain vs tail", sel ? `${sel.usable_pairs} gene pairs` : ""],
    [pMa && pMa.full.fold ? `${f(pMa.full.fold, 1)}×` : "–", "SOG1 motif in banana target promoters", pMa ? `vs random genes, p = ${pfmt(pMa.full.p_perm)}` : ""],
  ];
  $("#kpis").innerHTML = kpis.map(([v, l, s]) => `<div class="kpi"><div class="v">${esc(v ?? "–")}</div><div class="l">${esc(l)}</div><div class="arms small">${esc(s || "")}</div></div>`).join("");

  // ------------------------------------------------------------ family
  const order = species.map((s) => s.species);
  (function copies() {
    const W = 520, row = 22, left = 190, H = order.length * row + 30;
    const s = svg(W, H);
    const max = Math.max(3, ...order.map((sp) => clade.filter((c) => c.species === sp).length));
    const x = (v) => left + (v / max) * (W - left - 30);
    for (let i = 0; i <= max; i++) add(s, "line", { x1: x(i), x2: x(i), y1: 6, y2: H - 22, stroke: css("--grid") });
    for (let i = 0; i <= max; i++) add(s, "text", { x: x(i), y: H - 6, "text-anchor": "middle", "font-size": 11, class: "muted" }, i);
    order.forEach((sp, i) => {
      const y = 8 + i * row;
      const info = spBy[sp];
      add(s, "rect", { x: 4, y: y + 4, width: 9, height: 9, rx: 2, fill: linColor(info.lineage) });
      add(s, "text", { x: 18, y: y + 13, "font-size": 12, "font-weight": sp === D.meta.focus_species ? 700 : 400 }, info.label);
      let acc = 0;
      LABELS.forEach((l) => {
        const n = clade.filter((c) => c.species === sp && c.clade_label === l).length;
        if (!n) return;
        const r = add(s, "rect", { x: x(acc), y: y + 2, width: x(acc + n) - x(acc), height: row - 6, fill: labColor[l], rx: 2 });
        add(r, "title", {}, `${info.label}: ${n} ${l}-like`);
        acc += n;
      });
      add(s, "text", { x: x(acc) + 5, y: y + 14, "font-size": 11, class: "muted" }, acc || "0");
    });
    $("#c-copies").appendChild(s);
    $("#c-copies").insertAdjacentHTML("afterend", `<div class="legend">${LABELS.map((l) => `<span><i class="sw" style="background:${labColor[l]}"></i>${l}-like</span>`).join("")}</div>`);
  })();

  table($("#t-counts"), [{ t: "Genome" }, { t: "Genes", num: 1 }, { t: "NAC", num: 1 }, { t: "Published", num: 1 }, { t: "SOG1 clade", num: 1 }],
    order.map((sp) => { const c = counts[sp] || {}; return [esc(spBy[sp].label), (c.genes || 0).toLocaleString(), c.nac,
      c.published ? `${c.published} <span class="${c.in_range ? "ok" : "bad"}">${c.in_range ? "✓" : "✗"}</span>` : "", c.sog1_clade]; }));
  const noNac = order.filter((sp) => counts[sp] && counts[sp].sog1_clade === 0).map((sp) => spBy[sp].label);
  $("#family-note").textContent = `Every protein-coding gene in ${species.length} genomes was searched with the Pfam NAC (NAM) domain model. ` +
    `Banana has ${fc.nac} NAC genes, ${focus.length} of them in the SOG1 clade (${focusLab.map(([l, n]) => `${n} closest to ${l}`).join(", ")}). ` +
    (noNac.length ? `No SOG1-clade gene was found in ${noNac.join(", ")}, ` + "consistent with SOG1 being a land-plant innovation." : "");

  // ------------------------------------------------------------ tree
  function parseNewick(s) {
    let i = 0;
    function node() {
      const n = { children: [] };
      if (s[i] === "(") {
        i++;
        do { n.children.push(node()); } while (s[i] === "," && ++i);
        i++; // ')'
      }
      let lab = "";
      while (i < s.length && !",():;".includes(s[i])) lab += s[i++];
      n.label = lab;
      if (s[i] === ":") { i++; let b = ""; while (i < s.length && !",();".includes(s[i])) b += s[i++]; n.len = parseFloat(b); }
      return n;
    }
    return node();
  }
  if (D.tree && D.tree.newick) {
    const root = parseNewick(D.tree.newick.trim());
    const tips = [];
    (function walk(n, depth) {
      n.depth = depth + (n.len || 0);
      if (!n.children.length) tips.push(n);
      n.children.forEach((c) => walk(c, n.depth));
    })(root, 0);
    tips.forEach((t, i) => (t.y = i));
    (function sety(n) { if (n.children.length) { n.children.forEach(sety); n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2; } })(root);
    const maxD = Math.max(...tips.map((t) => t.depth));
    const W = 980, row = 17, H = tips.length * row + 40, left = 10, right = 330;
    const X = (d) => left + (d / maxD) * (W - left - right);
    const Y = (y) => 14 + y * row;
    const s = svg(W, H);
    const names = D.tree.names || {};
    (function draw(n) {
      n.children.forEach((c) => {
        add(s, "path", { d: `M${X(n.depth)},${Y(c.y)}H${X(c.depth)}`, stroke: css("--ink"), "stroke-width": 1.2, fill: "none" });
        draw(c);
      });
      if (n.children.length) {
        add(s, "path", { d: `M${X(n.depth)},${Y(n.children[0].y)}V${Y(n.children[n.children.length - 1].y)}`, stroke: css("--ink"), "stroke-width": 1.2 });
        const sup = parseFloat(n.label);
        if (isFinite(sup) && n !== root) add(s, "text", { x: X(n.depth) - 3, y: Y(n.y) - 3, "text-anchor": "end", "font-size": 9.5, class: "muted", "font-weight": sup >= 95 ? 600 : 400 }, Math.round(sup));
      } else {
        const key = names[n.label] || n.label.replace("__", "|");
        const c = byKey[key];
        const [sh, gene] = key.split("|");
        const sp = shortBy[sh] || {};
        const isOut = !c;
        const lab = `${sp.label ? sp.label.replace(/ \(.*\)/, "") : sh} · ${gene}` + (c && c.symbol ? ` (${c.symbol})` : isOut ? " (outgroup)" : "");
        add(s, "circle", { cx: X(n.depth), cy: Y(n.y), r: 3.4, fill: isOut ? css("--muted") : linColor(sp.lineage) });
        add(s, "text", { x: X(n.depth) + 7, y: Y(n.y) + 4, class: "tree-label" + (sh === focusShort ? " focus" : ""), fill: isOut ? css("--muted") : css("--ink") }, lab);
        if (c && c.clade_label) add(s, "text", { x: W - 8, y: Y(n.y) + 4, "text-anchor": "end", "font-size": 11, fill: labColor[c.clade_label] }, `${c.clade_label}-like`);
      }
    })(root);
    const sb = 0.2, sx = X(0), sw = X(sb) - X(0);
    add(s, "path", { d: `M${sx},${H - 14}h${sw}`, stroke: css("--ink") });
    add(s, "text", { x: sx + sw / 2, y: H - 2, "text-anchor": "middle", "font-size": 10, class: "muted" }, `${sb} substitutions/site`);
    $("#c-tree").appendChild(s);
    $("#tree-legend").innerHTML = [...new Set(species.map((s) => s.lineage))].map((l) => `<span><i class="sw" style="background:${linColor(l)}"></i>${l}</span>`).join("");
    $("#tree-note").textContent = `${D.tree.sequences} proteins, ${D.tree.trimmed_columns} of ${D.tree.aligned_columns} alignment columns kept after trimming; ` +
      `best model by BIC: ${D.tree.model}. The 2023 report used a neighbour-joining tree of five proteins; this tree uses every SOG1-clade gene found.`;
    window.__tipOrder = tips.map((t) => names[t.label] || t.label);
  } else {
    $("#tree-note").textContent = "The tree has not been built yet.";
  }

  // ------------------------------------------------------------ SQ maps
  const tipOrder = (window.__tipOrder || []).filter((k) => byKey[k]);
  const mapOrder = tipOrder.length ? tipOrder : clade.map((c) => c.key);
  (function maps() {
    const W = 900, row = 18, left = 210, H = mapOrder.length * row + 34;
    const maxL = Math.max(...mapOrder.map((k) => byKey[k].length));
    const x = (p) => left + (p / maxL) * (W - left - 20);
    const s = svg(W, H);
    for (let p = 0; p <= maxL; p += 100) {
      add(s, "line", { x1: x(p), x2: x(p), y1: 4, y2: H - 22, stroke: css("--grid") });
      add(s, "text", { x: x(p), y: H - 8, "text-anchor": "middle", "font-size": 10, class: "muted" }, p);
    }
    mapOrder.forEach((k, i) => {
      const c = byKey[k], y = 10 + i * row;
      const sp = shortBy[c.short];
      add(s, "text", { x: left - 8, y: y + 4, "text-anchor": "end", "font-size": 11, "font-weight": c.short === focusShort ? 700 : 400 }, `${c.short} ${c.gene}`);
      add(s, "line", { x1: x(0), x2: x(c.length), y1: y, y2: y, stroke: css("--muted"), "stroke-width": 2 });
      add(s, "rect", { x: x(c.domain_start), y: y - 5, width: x(c.domain_end) - x(c.domain_start), height: 10, rx: 3, fill: css("--accent-soft"), stroke: css("--accent") });
      arr(D.sq.positions[k]).forEach((p) => {
        const t = add(s, "line", { x1: x(p), x2: x(p), y1: y - 7, y2: y + 7, stroke: p > c.domain_end ? css("--c3") : css("--muted"), "stroke-width": 2 });
        add(t, "title", {}, `${c.short} ${c.gene}: ${byKey[k] && D.profiles[k] ? D.profiles[k].seq.substr(p - 1, 2) : "S/TQ"} at ${p}`);
      });
      if (c.scd) add(s, "text", { x: x(c.length) + 4, y: y + 4, "font-size": 10, fill: css("--c3") }, "SCD");
    });
    $("#c-maps").appendChild(s);
  })();
  const tv = D.sq.tail_vs_domain;
  $("#sq-note").textContent = `Each line is one SOG1-clade protein (same order as the tree): green box = NAC DNA-binding domain, pink ticks = SQ/TQ sites in the tail. ` +
    `Median density is ${f(tv.median_tail)} sites per 100 residues in the tail vs ${f(tv.median_domain)} in the domain (paired Wilcoxon p = ${pfmt(tv.wilcoxon_p)}, ${tv.n} proteins). ` +
    `${sqSig} of ${sqRows.length} tails have more sites than their own composition gives by chance (q < 0.05).`;
  table($("#t-sq"), [{ t: "Protein" }, { t: "Closest" }, { t: "Length", num: 1 }, { t: "Tail SQ/TQ", num: 1 }, { t: "Expected", num: 1 }, { t: "q", num: 1 }, { t: "SCD" }],
    mapOrder.map((k) => { const c = byKey[k]; return [`${c.short === focusShort ? "<b>" : ""}${esc(geneName(c))}${c.short === focusShort ? "</b>" : ""}`, c.clade_label,
      c.length, c.sq_tail, f(c.sq_tail_expected, 1), `<span class="${c.sq_tail_q < 0.05 ? "ok" : ""}">${pfmt(c.sq_tail_q)}</span>`, c.scd ? '<span class="tag">yes</span>' : '<span class="tag no">no</span>']; }));
  if (D.sq.conservation && D.sq.conservation.length) {
    table($("#t-cons"), [{ t: "Arabidopsis SOG1 site" }, { t: "In SOG1-like orthologs", num: 1 }].concat(focus.map((c) => ({ t: c.gene }))),
      D.sq.conservation.map((e) => [`${esc(e.motif)}${e.at_position}`, e.conserved_in_sog1_orthologs == null ? "–" : Math.round(100 * e.conserved_in_sog1_orthologs) + "%"]
        .concat(focus.map((c) => (e.present[c.key] ? '<span class="ok">✓</span>' : '<span class="bad">✗</span>')))));
  } else $("#cons-h").hidden = true;

  // ------------------------------------------------------------ disorder profiles
  const prof = Object.keys(D.profiles).filter((k) => byKey[k]);
  const pick = $("#dis-pick");
  pick.innerHTML = prof.map((k) => `<option value="${esc(k)}">${esc(geneName(byKey[k]))} (${byKey[k].clade_label}-like)</option>`).join("");
  const defaultKey = prof.find((k) => k.startsWith(focusShort + "|") && byKey[k].clade_label === "SOG1") || prof[0];
  if (defaultKey) pick.value = defaultKey;
  function drawDis() {
    const k = pick.value, c = byKey[k], p = D.profiles[k];
    const el = $("#c-dis"); el.innerHTML = "";
    if (!p) return;
    const W = 900, H = 230, L = 40, R = 10, T = 10, B = 30, n = p.seq.length;
    const x = (i) => L + (i / n) * (W - L - R), y = (v) => T + (1 - v) * (H - T - B);
    const s = svg(W, H);
    add(s, "rect", { x: x(c.domain_start - 1), y: T, width: x(c.domain_end) - x(c.domain_start - 1), height: H - T - B, fill: css("--accent-soft") });
    [0, 0.5, 1].forEach((v) => { add(s, "line", { x1: L, x2: W - R, y1: y(v), y2: y(v), stroke: css("--grid") }); add(s, "text", { x: L - 6, y: y(v) + 4, "text-anchor": "end", "font-size": 10, class: "muted" }, v); });
    for (let i = 0; i <= n; i += 100) add(s, "text", { x: x(i), y: H - 10, "text-anchor": "middle", "font-size": 10, class: "muted" }, i);
    const line = (vals, color) => add(s, "path", { d: vals.map((v, i) => `${i ? "L" : "M"}${x(i + 0.5).toFixed(1)},${y(v).toFixed(1)}`).join(""), stroke: color, fill: "none", "stroke-width": 1.6 });
    if (p.disorder) line(p.disorder, css("--c1"));
    if (p.plddt) line(p.plddt.map((v) => v / 100), css("--c2"));
    arr(D.sq.positions[k]).forEach((q) => add(s, "line", { x1: x(q - 0.5), x2: x(q - 0.5), y1: H - B, y2: H - B - 10, stroke: css("--c3"), "stroke-width": 2 }));
    el.appendChild(s);
  }
  pick.addEventListener("change", drawDis); drawDis();
  const af = D.sq.alphafold || [];
  const afAt = af.find((a) => a.key.startsWith("At|"));
  $("#dis-note").textContent = `SOG1 needs a folded NAC domain to bind DNA and a flexible tail for the kinases to reach. ` +
    (clade.some((c) => c.disorder_tail != null) ? `Across the clade, a median ${Math.round(100 * median(clade.map((c) => c.disorder_tail)))}% of tail residues are predicted disordered vs ${Math.round(100 * median(clade.map((c) => c.disorder_domain)))}% of domain residues. ` : "") +
    (af.length ? `AlphaFold models with an identical sequence were found for ${af.length} proteins; ${afAt && afAt.rho != null ? `for Arabidopsis SOG1 the disorder prediction and AlphaFold confidence agree (Spearman ρ = ${f(afAt.rho)}), ` : ""}` +
      `median ${Math.round(100 * median(af.map((a) => a.tail_low_plddt)))}% of tail residues have very low AlphaFold confidence (pLDDT < 50).` : "");
  function median(v) { v = v.filter((x) => x != null && isFinite(x)).sort((a, b) => a - b); return v.length ? (v.length % 2 ? v[(v.length - 1) / 2] : (v[v.length / 2 - 1] + v[v.length / 2]) / 2) : NaN; }

  // ------------------------------------------------------------ selection
  if (sel && D.pairs.length) {
    const W = 420, H = 380, L = 50, B = 44, T = 10, R = 10;
    const lo = Math.log10(0.01), hi = Math.log10(3);
    const tx = (v) => L + ((Math.log10(Math.max(v, 0.01)) - lo) / (hi - lo)) * (W - L - R);
    const ty = (v) => T + (1 - (Math.log10(Math.max(v, 0.01)) - lo) / (hi - lo)) * (H - T - B);
    const s = svg(W, H);
    [0.01, 0.1, 1].forEach((v) => {
      add(s, "line", { x1: tx(v), x2: tx(v), y1: T, y2: H - B, stroke: css("--grid") }); add(s, "line", { x1: L, x2: W - R, y1: ty(v), y2: ty(v), stroke: css("--grid") });
      add(s, "text", { x: tx(v), y: H - B + 16, "text-anchor": "middle", "font-size": 10, class: "muted" }, v);
      add(s, "text", { x: L - 6, y: ty(v) + 4, "text-anchor": "end", "font-size": 10, class: "muted" }, v);
    });
    add(s, "line", { x1: tx(0.01), y1: ty(0.01), x2: tx(3), y2: ty(3), stroke: css("--muted"), "stroke-dasharray": "4 4" });
    D.pairs.forEach((p) => {
      const ban = p.species_a === D.meta.focus_species || p.species_b === D.meta.focus_species;
      const c = add(s, "circle", { cx: tx(p.omega_domain), cy: ty(p.omega_tail), r: ban ? 3.6 : 2.6, fill: ban ? css("--c1") : css("--muted"), "fill-opacity": ban ? 0.85 : 0.45 });
      add(c, "title", {}, `${p.a} vs ${p.b}\nω domain ${f(p.omega_domain, 3)}, tail ${f(p.omega_tail, 3)}, dS ${f(p.ds_full)}`);
    });
    add(s, "text", { x: (L + W - R) / 2, y: H - 8, "text-anchor": "middle", "font-size": 12 }, "dN/dS of NAC domain");
    add(add(s, "text", { x: 14, y: (T + H - B) / 2, "text-anchor": "middle", "font-size": 12, transform: `rotate(-90 14 ${(T + H - B) / 2})` }), "tspan", {}, "dN/dS of tail");
    $("#c-sel").appendChild(s);
    $("#c-sel").insertAdjacentHTML("afterend", `<div class="legend"><span><i class="sw" style="background:${css("--c1")}"></i>pair with a banana gene</span><span><i class="sw" style="background:${css("--muted")}"></i>other pairs</span></div>`);
    $("#sel-note").textContent = `Across ${sel.usable_pairs} gene pairs with usable divergence, the median dN/dS is ${f(sel.median_omega_domain, 3)} in the NAC domain and ${f(sel.median_omega_tail, 3)} in the tail ` +
      `(tail higher in ${Math.round(100 * sel.tail_higher_fraction)}% of pairs; paired Wilcoxon p = ${pfmt(sel.wilcoxon_p)}). Both are far below 1: the whole protein is under purifying selection, the DNA-binding domain most strongly.`;
  } else $("#sel-note").textContent = "Selection results not available.";
  const para = (sel && sel.focus_paralogs) || [];
  table($("#t-para"), [{ t: "Pair" }, { t: "Closest" }, { t: "Ks", num: 1 }, { t: "dN/dS", num: 1 }],
    para.length ? para.map((p) => [`<span class="mono">${esc(p.a.split("|")[1])}<br>${esc(p.b.split("|")[1])}</span>`, `${p.label_a} / ${p.label_b}`, f(p.ks), f(p.omega, 3)]) : [["Only one banana gene in the clade", "", "", ""]]);

  // ------------------------------------------------------------ promoters
  (function prom() {
    const sps = Object.keys(P);
    const motifs = D.promoters.motifs;
    const W = 420, H = 260, L = 44, B = 50, T = 14;
    const groups = sps.flatMap((sp) => motifs.map((m) => ({ sp, m, r: P[sp][m] })));
    const max = Math.max(0.5, ...groups.map((g) => Math.max(g.r.full.mean_target || 0, g.r.full.mean_background || 0)));
    const y = (v) => T + (1 - v / max) * (H - T - B);
    const s = svg(W, H);
    [0, max / 2, max].forEach((v) => { add(s, "line", { x1: L, x2: W - 6, y1: y(v), y2: y(v), stroke: css("--grid") }); add(s, "text", { x: L - 6, y: y(v) + 4, "text-anchor": "end", "font-size": 10, class: "muted" }, f(v, 1)); });
    const gw = (W - L - 10) / sps.length;
    sps.forEach((sp, i) => {
      const x0 = L + i * gw + 10;
      motifs.forEach((m, j) => {
        const r = P[sp][m].full, bw = 22, xx = x0 + j * (2 * bw + 14);
        const ctl = j > 0;
        add(s, "rect", { x: xx, y: y(r.mean_target), width: bw, height: y(0) - y(r.mean_target), fill: ctl ? css("--muted") : css("--c1"), rx: 2 });
        add(s, "rect", { x: xx + bw, y: y(r.mean_background), width: bw, height: y(0) - y(r.mean_background), fill: ctl ? css("--line") : css("--accent-soft"), stroke: ctl ? css("--muted") : css("--accent"), rx: 2 });
        add(s, "text", { x: xx + bw, y: y(Math.max(r.mean_target, r.mean_background)) - 4, "text-anchor": "middle", "font-size": 10, class: ctl ? "muted" : "" }, `p ${pfmt(r.p_perm)}`);
      });
      add(s, "text", { x: x0 + 55, y: H - 26, "text-anchor": "middle", "font-size": 12, "font-weight": 600 }, spBy[sp].label.replace(/ \(.*\)/, ""));
      add(s, "text", { x: x0 + 55, y: H - 12, "text-anchor": "middle", "font-size": 10, class: "muted" }, `${P[sp][motifs[0]].full.n_target} targets · ${P[sp][motifs[0]].full.n_background} random`);
    });
    $("#c-prom").appendChild(s);
    $("#c-prom").insertAdjacentHTML("afterend", `<div class="legend"><span><i class="sw" style="background:${css("--c1")}"></i>targets</span><span><i class="sw" style="background:${css("--accent-soft")};border:1px solid ${css("--accent")}"></i>random genes</span></div>`);

    const seg = $("#dens-sp");
    seg.innerHTML = sps.map((sp, i) => `<button data-sp="${sp}" aria-pressed="${i === 0}">${esc(spBy[sp].label.replace(/ \(.*\)/, ""))}</button>`).join("");
    function dens(sp) {
      const el = $("#c-dens"); el.innerHTML = "";
      const r = P[sp][motifs[0]], n = r.density_target.length;
      const W2 = 420, H2 = 220, L2 = 40, B2 = 34, T2 = 10;
      const mx = Math.max(0.5, ...r.density_target, ...r.density_background);
      const x = (i) => L2 + (i / n) * (W2 - L2 - 8), yy = (v) => T2 + (1 - v / mx) * (H2 - T2 - B2);
      const s2 = svg(W2, H2);
      [0, mx / 2, mx].forEach((v) => { add(s2, "line", { x1: L2, x2: W2 - 8, y1: yy(v), y2: yy(v), stroke: css("--grid") }); add(s2, "text", { x: L2 - 6, y: yy(v) + 4, "text-anchor": "end", "font-size": 10, class: "muted" }, f(v, 1)); });
      const bw = (W2 - L2 - 8) / n;
      r.density_background.forEach((v, i) => add(s2, "rect", { x: x(i) + 1, y: yy(v), width: bw - 2, height: yy(0) - yy(v), fill: css("--accent-soft"), stroke: css("--accent") }));
      add(s2, "path", { d: r.density_target.map((v, i) => `${i ? "L" : "M"}${x(i + 0.5)},${yy(v)}`).join(""), stroke: css("--c1"), "stroke-width": 2, fill: "none" });
      r.density_target.forEach((v, i) => add(s2, "circle", { cx: x(i + 0.5), cy: yy(v), r: 2.5, fill: css("--c1") }));
      for (let i = 0; i <= n; i += 5) add(s2, "text", { x: x(i), y: H2 - 18, "text-anchor": "middle", "font-size": 10, class: "muted" }, `${i * 100}`);
      add(s2, "text", { x: (L2 + W2) / 2, y: H2 - 3, "text-anchor": "middle", "font-size": 11, class: "muted" }, "bp upstream of gene start");
      el.appendChild(s2);
    }
    seg.addEventListener("click", (e) => { const b = e.target.closest("button"); if (!b) return; seg.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", x === b)); dens(b.dataset.sp); });
    if (sps.length) dens(sps[0]);

    const at = P.arabidopsis_thaliana, fo = P[D.meta.focus_species];
    const atR = at && at[motifs[0]], foR = fo && fo[motifs[0]];
    $("#prom-note").textContent = (atR ? `Positive control: promoters of ${atR.full.n_target} reported SOG1 targets in Arabidopsis carry on average ${f(atR.full.mean_target)} CTT(N)7AAG sites vs ${f(atR.full.mean_background)} in random genes (${f(atR.full.fold, 1)}×, p = ${pfmt(atR.full.p_perm)}). ` : "") +
      (foR ? `In banana, the ${foR.full.n_target} orthologs of those targets carry ${f(foR.full.mean_target)} vs ${f(foR.full.mean_background)} (${f(foR.full.fold, 1)}×, p = ${pfmt(foR.full.p_perm)}). ` : "") +
      `The control motif with the same letters tells whether a difference is about the SOG1 site or just base composition.`;

    const tg = D.targets;
    const atGenes = Object.fromEntries((at ? at.genes : []).map((g) => [g.gene, g]));
    const foGenes = Object.fromEntries((fo ? fo.genes : []).map((g) => [g.gene, g]));
    const rows = [];
    const seen = new Set();
    tg.forEach((t) => {
      const a = atGenes[t.at_gene];
      const firstRow = !seen.has(t.at_gene); seen.add(t.at_gene);
      const b = t.focus_gene ? foGenes[t.focus_gene] : null;
      const pos = (g) => (g ? (g.positions.length ? g.positions.slice(0, 6).map((p) => `−${p}`).join(", ") : "none") : "–");
      rows.push([firstRow ? `<b>${esc(t.symbol)}</b>` : "", firstRow ? `<span class="mono">${t.at_gene}</span>` : "", firstRow ? pos(a) : "",
        t.focus_gene ? `<span class="mono">${esc(t.focus_gene)}</span>` : '<span class="small">no reciprocal hit</span>', b ? pos(b) : "–"]);
    });
    table($("#t-targets"), [{ t: "Gene" }, { t: "Arabidopsis" }, { t: "Motif sites" }, { t: "Banana ortholog" }, { t: "Motif sites" }], rows);
  })();

  // ------------------------------------------------------------ checks
  const checks = [];
  Object.entries(D.counts).forEach(() => {});
  D.counts.filter((c) => c.published).forEach((c) => checks.push([c.in_range, `${c.label}: ${c.nac} NAC genes found, ${c.published} published`]));
  const atClade = clade.filter((c) => c.short === "At").map((c) => c.gene).sort();
  checks.push([atClade.join() === Object.keys(D.meta.at_clade).sort().join(), `Arabidopsis SOG1, NAC044 and NAC085 recovered, and nothing else (${atClade.length} genes)`]);
  const cr = counts.chlamydomonas_reinhardtii;
  if (cr) checks.push([cr.sog1_clade === 0, `No SOG1-clade gene in the green alga Chlamydomonas (negative control): ${cr.sog1_clade}`]);
  if (D.compara) {
    const both = D.compara.reduce((a, r) => a + r.both, 0), hmm = D.compara.reduce((a, r) => a + r.hmm_only.length, 0), cmp = D.compara.reduce((a, r) => a + r.compara_only.length, 0);
    checks.push([both / Math.max(1, both + hmm + cmp) >= 0.8, `Clade agrees with Ensembl Compara for ${both} genes; ${hmm} found only by my search, ${cmp} only by Compara`]);
  }
  if (sel) checks.push([sel.cds_mismatch.length <= 0.1 * (sel.codon_alignments + sel.cds_mismatch.length), `${sel.codon_alignments} coding sequences translate exactly to their protein (${sel.cds_mismatch.length} did not and were left out of dN/dS)`]);
  if (afAt && afAt.rho != null) checks.push([afAt.rho < -0.3, `Disorder prediction agrees with AlphaFold for Arabidopsis SOG1 (ρ = ${f(afAt.rho)})`]);
  if (atR_ok()) checks.push([P.arabidopsis_thaliana[M0].full.p_perm < 0.05, `Positive control: SOG1 motif enriched in Arabidopsis target promoters (p = ${pfmt(P.arabidopsis_thaliana[M0].full.p_perm)})`]);
  function atR_ok() { return P.arabidopsis_thaliana && P.arabidopsis_thaliana[M0]; }
  $("#checklist").innerHTML = checks.map(([ok, t]) => `<li class="${ok ? "ok" : "bad"}">${esc(t)}</li>`).join("");
  if (D.compara) table($("#t-compara"), [{ t: "Genome" }, { t: "Both", num: 1 }, { t: "Only my search" }, { t: "Only Compara" }],
    D.compara.map((r) => [esc(spBy[r.species] ? spBy[r.species].label : r.species), r.both, `<span class="mono">${r.hmm_only.join(", ") || "–"}</span>`, `<span class="mono">${r.compara_only.join(", ") || "–"}</span>`]));
  else $("#t-compara").outerHTML = '<p class="small">Ensembl Compara was not reachable during this run.</p>';

  // ------------------------------------------------------------ takeaways and phase 3
  const take = [];
  take.push(`Banana carries ${focus.length} SOG1-clade genes (${focusLab.map(([l, n]) => `${n} ${l}-like`).join(", ")}), against ${counts.arabidopsis_thaliana ? counts.arabidopsis_thaliana.sog1_clade : 3} in Arabidopsis and ${counts.oryza_sativa ? counts.oryza_sativa.sog1_clade : "?"} in rice.`);
  if (sqRows.length) take.push(`The ATM/ATR switch is ancient: ${sqSig} of ${sqRows.length} clade tails hold more SQ/TQ sites than chance and ${scd} carry an SQ/TQ cluster, banana included.`);
  if (sel && sel.median_omega_domain != null) take.push(`The DNA-binding domain is the most constrained part (median dN/dS ${f(sel.median_omega_domain, 3)} vs ${f(sel.median_omega_tail, 3)} in the tail).`);
  if (atR_ok() && pMa) take.push(`The SOG1 binding site marks target promoters in Arabidopsis (${f(P.arabidopsis_thaliana[M0].full.fold, 1)}×) ${pMa.full.p_perm < 0.05 ? "and also in their banana orthologs" : "but not clearly in their banana orthologs"} (${f(pMa.full.fold, 1)}×, p = ${pfmt(pMa.full.p_perm)}).`);
  $("#takeaways").innerHTML = take.map((t) => `<p>${esc(t)}</p>`).join("");

  const fg = D.meta.focus_gene;
  $("#next-list").innerHTML = [
    `Complementation: transform the Arabidopsis <i>sog1-1</i> mutant with the pRI101-AN construct from phase 1${fg ? ` (<span class="mono">${esc(fg)}</span>)` : ""} and test root growth and cell death after gamma rays or zeocin.`,
    `Phospho-dead control: mutate the conserved tail SQ sites to AQ and check that the rescue is lost, as shown for Arabidopsis SOG1.`,
    `Paralogs: measure each banana SOG1-clade copy by qPCR after DNA damage to see which copy responds.`,
    `Targets: ChIP-qPCR or a luciferase reporter on banana target promoters with the strongest CTT(N)7AAG sites.`,
  ].map((t) => `<li>${t}</li>`).join("");
})();
