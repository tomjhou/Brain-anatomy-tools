"""The small script every 3D page from the ABC browser carries so that injection sites (and optic fibers) can be shown in it without the page being
rewritten: the page is opened with  #sites=<address of a side file>  after its own address, loads that side file (a little .js file written by
add_injections_3d.py that calls window.__injectionSites({...})) and draws what is in it. Opened without the #sites=, the page shows no sites. A tab keeps
what it was opened with: a refresh loads the same side file, and a later export opens a new tab with its own side file, so the sets can be compared.

The script is used in two places: the ABC browser puts it in every page it writes (write_3d_group_html), and compressed_page.upgrade_page adds it to
pages written before that."""

LOADER_MARK = 'injection-sites-loader'

LOADER_JS = r"""
/* injection-sites-loader: draws the injection sites of the side file named after # in the address (#sites=...), if any */
(function () {
  var match = /[#&]sites=([^&]*)/.exec(location.hash);
  if (!match) return;
  var source = decodeURIComponent(match[1]);
  // only a file on this computer (file://...) or one beside the page (name.js): a link that sends a page somewhere else gets nothing
  if (!/^(file:\/\/|[A-Za-z0-9_.-]+\.js$)/.test(source)) return;
  function note(text) {
    var box = document.createElement('div');
    box.textContent = text;
    box.style.cssText = 'position:fixed;left:10px;bottom:10px;z-index:30;background:#202530;color:#e6e8ee;border:1px solid #5a6178;border-radius:4px;' +
                        'padding:6px 10px;font:12px sans-serif;max-width:60%;';
    document.body.appendChild(box);
    setTimeout(function () { box.remove(); }, 12000);
  }
  function draw(data) {
    var gd = document.querySelector('.plotly-graph-div');
    if (!window.Plotly || !gd || !gd._fullLayout || !gd._fullLayout.scene) return setTimeout(function () { draw(data); }, 200);
    if (data.label) document.title = 'Injections: ' + data.label;
    var uv = data.sphere.v, uf = data.sphere.f, n = uv.length / 3;
    var traces = [];
    data.groups.forEach(function (g) {
      var sites = data.sites.filter(function (s) { return s[6] === g.kind; });
      if (!sites.length) return;
      var r = g.radius, x = [], y = [], z = [], i = [], j = [], k = [], text = [];
      sites.forEach(function (s, b) {
        for (var v = 0; v < n; v++) {
          x.push(s[0] + r * uv[3 * v]); y.push(s[1] + r * uv[3 * v + 1]); z.push(s[2] + r * uv[3 * v + 2]);
          text.push(s[3]);
        }
        for (var f = 0; f < uf.length; f += 3) { i.push(b * n + uf[f]); j.push(b * n + uf[f + 1]); k.push(b * n + uf[f + 2]); }
      });
      // the tracks: from where the DV was zeroed (brain surface, or bregma's height) down to the ball
      var tx = [], ty = [], tz = [], ttext = [];
      sites.forEach(function (s) {
        tx.push(s[4], s[0], null); ty.push(s[1], s[1], null); tz.push(s[5], s[2], null); ttext.push(s[3], s[3], null);
      });
      traces.push(
        {type: 'mesh3d', x: x, y: y, z: z, i: i, j: j, k: k, color: g.color, opacity: 1,
         flatshading: false, lighting: {ambient: 0.55, diffuse: 0.8, specular: 0.3, roughness: 0.5},
         name: g.name + ' (' + sites.length + ')', showlegend: true, text: text, hoverinfo: 'text'},
        {type: 'scatter3d', mode: 'lines', x: tx, y: ty, z: tz, line: {color: g.color, width: g.width},
         connectgaps: false, name: g.trackName, showlegend: true, text: ttext, hoverinfo: 'text'});
    });
    Plotly.addTraces(gd, traces);
    if (!gd.__hoverPassThrough) passThrough(gd);   // (pages from before the ABC browser did this itself)
  }
  // Hover through the background, for pages written before the ABC browser did it (the same code as its THREED_WEB_LEGEND_ROUTER_SCRIPT): objects whose
  // trace skips hover are left out of plotly's pick pass, so the translucent brain surface no longer hides the balls and cells inside it from the pointer.
  function passThrough(gd) {
    gd.__hoverPassThrough = function () {
      var scene = gd._fullLayout && gd._fullLayout.scene && gd._fullLayout.scene._scene;
      if (!scene || !scene.traces || !scene.glplot) return;
      var changed = false;
      Object.keys(scene.traces).forEach(function (uid) {
        var trace = scene.traces[uid];
        if (!trace || !trace.data || trace.data.hoverinfo !== 'skip') return;
        Object.keys(trace).forEach(function (key) {
          var object = trace[key];
          if (object && typeof object === 'object' && typeof object.drawPick === 'function') { object.drawPick = null; changed = true; }
        });
      });
      if (changed) scene.glplot.update();
    };
    ['plotly_afterplot', 'plotly_restyle'].forEach(function (name) {
      gd.on(name, function () { setTimeout(gd.__hoverPassThrough, 0); });
    });
    gd.__hoverPassThrough();
  }
  window.__injectionSites = draw;
  var tag = document.createElement('script');
  tag.src = source;
  tag.onerror = function () { note('The injection sites file could not be loaded (it may have been cleaned up): ' + source); };
  document.head.appendChild(tag);
})();
"""
