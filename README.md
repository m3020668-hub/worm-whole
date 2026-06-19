# worm-whole

<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Wormhole Simulator — Earth ↔ Andromeda (Morris-Thorne GR)</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<style>
* { margin:0; padding:0; box-sizing:border-box; }
html,body { width:100%; height:100%; overflow:hidden; background:#000308; }
canvas { display:block; }

#ui {
  position:fixed; top:14px; left:14px; width:268px;
  background:rgba(3,8,22,0.90); border:1px solid #1a3660;
  border-radius:10px; padding:13px 14px;
  font-family:'Courier New',monospace; color:#b8d4ee; font-size:11px;
  z-index:20; user-select:none; backdrop-filter:blur(6px);
}
#ui h2 { font-size:12px; letter-spacing:1.5px; color:#4db6ff;
          margin-bottom:10px; text-transform:uppercase; }
.eq-box { background:rgba(0,15,45,0.55); border-left:2px solid #1e4080;
           padding:5px 8px; margin-bottom:9px;
           font-size:10px; color:#5585a8; line-height:1.8; }
.sep { border:none; border-top:1px solid #162840; margin:9px 0; }
label { display:flex; justify-content:space-between; align-items:center;
        margin:5px 0; font-size:10.5px; color:#8ab0cc; }
label span.n { flex:1; }
input[type=range] { width:98px; accent-color:#2a6ad5; cursor:pointer; }
.v { color:#4db6ff; min-width:40px; text-align:right; font-size:10.5px; }
.brow { display:flex; flex-wrap:wrap; gap:4px; margin:3px 0; }
button {
  flex:1 1 45%; padding:4px 5px; font-size:10px;
  background:rgba(16,38,90,0.65); color:#88b0d0;
  border:1px solid #1a3660; border-radius:5px;
  cursor:pointer; font-family:'Courier New',monospace; transition:background .15s;
}
button:hover { background:rgba(35,75,180,0.75); color:#fff; }
button.off  { background:rgba(6,12,28,0.6); color:#304050; border-color:#0d1e30; }
button.act  { background:rgba(45,105,220,0.75); color:#ddeeff; border-color:#3a78ee; }
.sl { display:flex; justify-content:space-between; margin:3px 0; font-size:10px; }
.sl .lbl { color:#44606e; }
.sl .val { color:#4db6ff; }
.neg { color:#ff5555 !important; }
.pos { color:#55ff88 !important; }

#panel-right {
  position:fixed; top:14px; right:14px; width:208px;
  background:rgba(3,8,22,0.85); border:1px solid #1a3660;
  border-radius:10px; padding:12px 13px;
  font-family:'Courier New',monospace; font-size:10px;
  color:#4a6070; line-height:1.75; z-index:20;
}
#panel-right h3 { font-size:11px; color:#4db6ff; margin-bottom:7px; letter-spacing:1px; }
.exotic  { color:#ff6050; }
.embed   { color:#50d4a0; }
.crit    { color:#88dd44; }
.legend  { margin-top:8px; }
.legend-row { display:flex; align-items:center; gap:6px; margin:3px 0; font-size:9.5px; color:#486070; }
.dot { width:10px; height:10px; border-radius:50%; flex-shrink:0; }

.hint {
  position:fixed; bottom:14px; left:50%; transform:translateX(-50%);
  background:rgba(3,8,22,0.72); border:1px solid #1a3660; border-radius:20px;
  padding:5px 18px; font-family:'Courier New',monospace;
  font-size:10px; color:#384858; z-index:20; white-space:nowrap;
}
</style>
</head>
<body>
<canvas id="c"></canvas>

<!-- LEFT PANEL -->
<div id="ui">
  <h2>⚛ Wormhole Simulator</h2>
  <div class="eq-box">
    ds² = −dt² + dl² + r(l)² dΩ²<br>
    r(l) = √(l² + b₀²) <em style="color:#33507a">(Morris-Thorne 1988)</em>
  </div>

  <label><span class="n">Throat radius b₀</span>
    <input type="range" id="sB0" min="0.5" max="3.5" step="0.1" value="1.5">
    <span class="v" id="vB0">1.5</span>
  </label>
  <label><span class="n">Tunnel length</span>
    <input type="range" id="sLen" min="3" max="10" step="0.5" value="6">
    <span class="v" id="vLen">6×b₀</span>
  </label>
  <label><span class="n">Anim speed</span>
    <input type="range" id="sSpd" min="0" max="3" step="0.1" value="1">
    <span class="v" id="vSpd">1.0×</span>
  </label>

  <hr class="sep">
  <div class="brow">
    <button id="btnWire">Wireframe OFF</button>
    <button id="btnRays">Light Rays ON</button>
    <button id="btnTravel">Traveler ON</button>
    <button id="btnGalaxy">Galaxies ON</button>
  </div>

  <hr class="sep">
  <div class="sl"><span class="lbl">ρ(throat) [geom. units]</span></div>
  <div class="sl"><span class="val neg" id="svRho">—</span></div>
  <div class="sl">
    <span class="lbl">NEC: ρ + p_r at throat</span>
    <span class="val neg">negative ✗</span>
  </div>
  <div class="sl">
    <span class="lbl">Tidal-safe crossing speed</span>
  </div>
  <div class="sl"><span class="val pos" id="svVmax">—</span></div>
</div>

<!-- RIGHT PANEL -->
<div id="panel-right">
  <h3>GR FIELD EQUATIONS</h3>
  G<sub>μν</sub> = 8π T<sub>μν</sub> demands:<br>
  <b class="exotic">ρ(l) = −b₀² / [8π(l²+b₀²)²]</b><br>
  Always negative → exotic matter. The Null Energy Condition is violated.<br><br>
  Embedding height (exact):<br>
  <b class="embed">z(l) = b₀ · arcsinh(l/b₀)</b><br><br>
  Critical ray impact param:<br>
  <b class="crit">b<sub>crit</sub> = b₀ · E</b><br>
  b &lt; b<sub>crit</sub> → transmitted ✓<br>
  b &gt; b<sub>crit</sub> → reflected ✗

  <div class="legend">
    <div class="legend-row"><div class="dot" style="background:#ffd700"></div>b=0 axial (transmitted)</div>
    <div class="legend-row"><div class="dot" style="background:#ff9900"></div>b=0.5b₀ (transmitted)</div>
    <div class="legend-row"><div class="dot" style="background:#ff4400"></div>b=0.88b₀ near-crit (transmitted)</div>
    <div class="legend-row"><div class="dot" style="background:#44ff88"></div>b=1.08b₀ near-crit (reflected)</div>
    <div class="legend-row"><div class="dot" style="background:#cc44ff"></div>b=1.6b₀ (reflected)</div>
    <div class="legend-row"><div class="dot" style="background:#00eeff"></div>traveler particle</div>
  </div>
</div>

<div class="hint">🖱 Drag: orbit &nbsp;|&nbsp; Right-drag / two-finger: pan &nbsp;|&nbsp; Scroll / pinch: zoom</div>

<script>
// ════════════════════════════════════════════════════════════════════
//  RENDERER + SCENE + CAMERA
// ════════════════════════════════════════════════════════════════════
const canvas = document.getElementById('c');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0x000410);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.01, 2000);

// Spherical-coordinate orbit state
let camTheta = -0.55, camPhi = 1.05, camRadius = 22;
let panX = 0, panY = 0;

function updateCamera() {
  const sx = Math.sin(camPhi), cx = Math.cos(camPhi);
  const st = Math.sin(camTheta), ct = Math.cos(camTheta);
  camera.position.set(
    panX + camRadius * sx * ct,
    panY + camRadius * cx,
    camRadius * sx * st
  );
  camera.lookAt(panX, panY, 0);
}
updateCamera();

// ════════════════════════════════════════════════════════════════════
//  LIGHTING
// ════════════════════════════════════════════════════════════════════
scene.add(new THREE.AmbientLight(0x08112a, 3.2));
const sun = new THREE.DirectionalLight(0x7799bb, 1.4);
sun.position.set(8, 14, 6);
scene.add(sun);
const throatGlow = new THREE.PointLight(0x2244ff, 4.0, 18);
scene.add(throatGlow);  // stays at origin (the throat)

// ════════════════════════════════════════════════════════════════════
//  PHYSICS PARAMETERS  (mutable via sliders)
// ════════════════════════════════════════════════════════════════════
let B0 = 1.5;    // throat radius (scene units; 1 unit ≈ 10 m for stat display)
let L_MAX = 6.0; // half-length of displayed tunnel in units of B0

// GR helper functions — exact closed forms from the metric
function rL(l, b)  { return Math.sqrt(l*l + b*b); }
function zL(l, b)  { return b * Math.asinh(l / b); }

// ════════════════════════════════════════════════════════════════════
//  WORMHOLE SURFACE  (THREE.LatheGeometry from the embedding)
//  LatheGeometry revolves a 2-D profile (x=radius, y=height) around Y.
//  For our metric: x = r(l) = √(l²+b₀²),  y = z(l) = b₀·arcsinh(l/b₀)
// ════════════════════════════════════════════════════════════════════
let wormMesh = null, wireMesh = null, wireMode = false;

function buildWormhole() {
  if (wormMesh) { scene.remove(wormMesh); wormMesh.geometry.dispose(); wormMesh.material.dispose(); }
  if (wireMesh) { scene.remove(wireMesh); wireMesh.geometry.dispose(); wireMesh.material.dispose(); }

  const N = 140, phiSegs = 90;
  const lMax = L_MAX * B0;
  const profile = [];
  for (let i = 0; i <= N; i++) {
    const l = -lMax + 2 * lMax * i / N;
    profile.push(new THREE.Vector2(rL(l, B0), zL(l, B0)));
  }

  const geo = new THREE.LatheGeometry(profile, phiSegs);

  // Solid surface
  wormMesh = new THREE.Mesh(geo, new THREE.MeshPhongMaterial({
    color: 0x143888, emissive: 0x05102a, specular: 0x3366cc,
    shininess: 65, transparent: true,
    opacity: wireMode ? 0.10 : 0.66, side: THREE.DoubleSide
  }));
  scene.add(wormMesh);

  // Wireframe overlay
  wireMesh = new THREE.Mesh(geo.clone(), new THREE.MeshBasicMaterial({
    color: 0x1e55bb, wireframe: true, transparent: true, opacity: 0.16
  }));
  wireMesh.visible = wireMode;
  scene.add(wireMesh);
}

// ════════════════════════════════════════════════════════════════════
//  THROAT RING  (red torus at l=0, where y=z_embed(0)=0)
// ════════════════════════════════════════════════════════════════════
let throatRingMesh = null;
function buildThroatRing() {
  if (throatRingMesh) { scene.remove(throatRingMesh); throatRingMesh.geometry.dispose(); }
  const geo = new THREE.TorusGeometry(B0, B0 * 0.022, 14, 120);
  throatRingMesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: 0xff2020 }));
  throatRingMesh.rotation.x = Math.PI / 2;  // lie in the XZ plane at y=0
  scene.add(throatRingMesh);
}

// ════════════════════════════════════════════════════════════════════
//  STARS  (built once; vertex-colored for realism)
// ════════════════════════════════════════════════════════════════════
(function buildStars() {
  const N = 7000;
  const pos = new Float32Array(N * 3), col = new Float32Array(N * 3);
  const palettes = [[1,.95,.8],[.8,.9,1],[1,1,1],[.9,.8,1]];
  for (let i = 0; i < N; i++) {
    const r = 450 + Math.random() * 550;
    const t = Math.random() * Math.PI * 2, p = Math.acos(2*Math.random()-1);
    pos[i*3]   = r*Math.sin(p)*Math.cos(t);
    pos[i*3+1] = r*Math.cos(p);
    pos[i*3+2] = r*Math.sin(p)*Math.sin(t);
    const c = palettes[Math.floor(Math.random()*4)];
    const b = 0.4 + Math.random()*0.6;
    col[i*3]=c[0]*b; col[i*3+1]=c[1]*b; col[i*3+2]=c[2]*b;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  geo.setAttribute('color',    new THREE.BufferAttribute(col, 3));
  scene.add(new THREE.Points(geo,
    new THREE.PointsMaterial({ size:0.65, sizeAttenuation:true, vertexColors:true })));
})();

// ════════════════════════════════════════════════════════════════════
//  LIGHT-RAY GEODESICS
//
//  Null-geodesic equations of motion in the equatorial plane
//  (derived from the metric by variation of the affine action):
//    dl/dλ   = v_l
//    dv_l/dλ = b² · l / r(l)⁴         (b = impact parameter Λ)
//    dφ/dλ   = b / r(l)²
//
//  Conserved quantities E and Λ set the initial v_l via
//    v_l₀ = √(E² − Λ²/r₀²)
//
//  Critical parameter b_crit = E·b₀:
//    b < b_crit  →  ray crosses the throat (transmitted)
//    b > b_crit  →  turning point before the throat (reflected)
// ════════════════════════════════════════════════════════════════════
function integrateRay(b, E, l0, steps, h, b0) {
  const r0 = rL(l0, b0);
  const disc = E*E - (b*b)/(r0*r0);
  if (disc <= 1e-14) return null;

  let l = l0, vl = Math.sqrt(disc), phi = 0;
  const path = [[l, phi]];

  for (let i = 0; i < steps; i++) {
    // RK4 integration of the three ODEs
    function D(ll, vvl) {
      const r2 = ll*ll + b0*b0, r4 = r2*r2;
      return [vvl, b*b*ll/r4, b/r2];
    }
    const [a1,b1,c1] = D(l, vl);
    const [a2,b2,c2] = D(l+h*.5*a1, vl+h*.5*b1);
    const [a3,b3,c3] = D(l+h*.5*a2, vl+h*.5*b2);
    const [a4,b4,c4] = D(l+h*a3,    vl+h*b3);
    l   += h/6*(a1+2*a2+2*a3+a4);
    vl  += h/6*(b1+2*b2+2*b3+b4);
    phi += h/6*(c1+2*c2+2*c3+c4);
    path.push([l, phi]);

    if (vl < 0 && l < l0*0.55) break;      // reflected ray headed back
    if (l >  Math.abs(l0)*1.7)   break;    // transmitted ray far enough
  }
  return path;
}

function pathTo3D(path, b0) {
  // (l, φ) → 3-D embedding: x=r·cos φ, y=z_embed, z=r·sin φ
  // This matches LatheGeometry's axis convention (revolves around Y).
  return path.map(([l, phi]) =>
    new THREE.Vector3(rL(l,b0)*Math.cos(phi), zL(l,b0), rL(l,b0)*Math.sin(phi)));
}

const RAY_DEFS = [
  { ratio:0.00, color:0xffd700 },   // axial (b=0), always transmitted
  { ratio:0.50, color:0xff9900 },   // transmitted
  { ratio:0.88, color:0xff4400 },   // near-critical, transmitted
  { ratio:1.08, color:0x44ff88 },   // near-critical, reflected
  { ratio:1.60, color:0xcc44ff },   // reflected
];

let rayLines = [], raysOn = true;

function buildRays() {
  for (const ln of rayLines) { scene.remove(ln); ln.geometry.dispose(); }
  rayLines = [];

  const l0 = -L_MAX * B0 * 0.97;
  const h  = Math.abs(l0) * 2.3 / 2800;

  for (const def of RAY_DEFS) {
    const path = integrateRay(def.ratio * B0, 1.0, l0, 2800, h, B0);
    if (!path) continue;
    const geo = new THREE.BufferGeometry().setFromPoints(pathTo3D(path, B0));
    const ln  = new THREE.Line(geo,
      new THREE.LineBasicMaterial({ color:def.color, transparent:true, opacity:0.82 }));
    ln.visible = raysOn;
    scene.add(ln);
    rayLines.push(ln);
  }
}

// ════════════════════════════════════════════════════════════════════
//  TRAVELER PARTICLE  (cyan sphere + animated trail)
// ════════════════════════════════════════════════════════════════════
const TRAIL = 110;
let travelerMesh = null, travelerTrail = null;
let trailBuf = null, trailAttr = null;
let travelT = 0, travelerOn = true;

function buildTraveler() {
  travelT = 0;
  if (travelerMesh)  { scene.remove(travelerMesh);  travelerMesh.geometry.dispose(); }
  if (travelerTrail) { scene.remove(travelerTrail); travelerTrail.geometry.dispose(); }

  travelerMesh = new THREE.Mesh(
    new THREE.SphereGeometry(B0 * 0.10, 14, 10),
    new THREE.MeshPhongMaterial({ color:0x00eeff, emissive:0x006688, shininess:90 })
  );
  travelerMesh.visible = travelerOn;
  scene.add(travelerMesh);

  trailBuf  = new Float32Array(TRAIL * 3);
  const startL = -L_MAX * B0;
  const sx = rL(startL, B0), sy = zL(startL, B0);
  for (let i = 0; i < TRAIL; i++) { trailBuf[i*3]=sx; trailBuf[i*3+1]=sy; trailBuf[i*3+2]=0; }
  trailAttr = new THREE.BufferAttribute(trailBuf, 3);
  const trailGeo = new THREE.BufferGeometry();
  trailGeo.setAttribute('position', trailAttr);
  travelerTrail = new THREE.Line(trailGeo,
    new THREE.LineBasicMaterial({ color:0x00cccc, transparent:true, opacity:0.50 }));
  travelerTrail.visible = travelerOn;
  scene.add(travelerTrail);
}

function stepTraveler(dt, speed) {
  travelT += dt * speed * 0.28;
  const lMax  = L_MAX * B0;
  const cycle = lMax * 2;
  const l     = -lMax + (((travelT % cycle) + cycle) % cycle);
  const px = rL(l, B0), py = zL(l, B0), pz = 0;  // travels at φ=0
  travelerMesh.position.set(px, py, pz);
  travelerMesh.scale.setScalar(0.85 + 0.22*Math.sin(travelT*3.4));

  // Shift trail ring-buffer
  for (let i = TRAIL-1; i > 0; i--) {
    trailBuf[i*3]   = trailBuf[(i-1)*3];
    trailBuf[i*3+1] = trailBuf[(i-1)*3+1];
    trailBuf[i*3+2] = trailBuf[(i-1)*3+2];
  }
  trailBuf[0]=px; trailBuf[1]=py; trailBuf[2]=pz;
  trailAttr.needsUpdate = true;
}

// ════════════════════════════════════════════════════════════════════
//  GALAXY HALOS  (disc of points at each wormhole mouth)
// ════════════════════════════════════════════════════════════════════
let galaxyObjs = [], galaxiesOn = true;

const lblEarth = document.createElement('div');
lblEarth.style.cssText = 'position:fixed;font:bold 11px Courier New;color:#88ff88;text-shadow:0 0 8px #00ee00;pointer-events:none;z-index:30;';
lblEarth.innerHTML = '🌍 MILKY WAY / EARTH';
document.body.appendChild(lblEarth);

const lblAndro = document.createElement('div');
lblAndro.style.cssText = 'position:fixed;font:bold 11px Courier New;color:#ffaadd;text-shadow:0 0 8px #ff44aa;pointer-events:none;z-index:30;';
lblAndro.innerHTML = '🌌 ANDROMEDA (M31)<br><span style="font-size:9px;color:#aa6688;font-weight:normal">2.537 million light-years</span>';
document.body.appendChild(lblAndro);

function buildGalaxies() {
  for (const g of galaxyObjs) {
    scene.remove(g);
    if (g.geometry) g.geometry.dispose();
    if (g.material) g.material.dispose();
  }
  galaxyObjs = [];

  const lMax = L_MAX * B0;
  const yTop = zL(lMax, B0), yBot = zL(-lMax, B0);

  function disc(cy, color, N, spread) {
    const pos = new Float32Array(N*3);
    for (let i = 0; i < N; i++) {
      const ang = Math.random()*Math.PI*2;
      const rr  = Math.pow(Math.random(), 0.38) * spread;
      pos[i*3]   = rr*Math.cos(ang);
      pos[i*3+1] = (Math.random()-.5)*spread*0.07;
      pos[i*3+2] = rr*Math.sin(ang);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    const pts = new THREE.Points(geo,
      new THREE.PointsMaterial({ color, size:0.26, sizeAttenuation:true, transparent:true, opacity:0.72 }));
    pts.position.set(0, cy, 0);
    pts.rotation.z = (Math.random()-.5)*0.55;
    pts.visible = galaxiesOn;
    scene.add(pts);
    galaxyObjs.push(pts);
  }

  disc(yBot - B0*5, 0x88bbff, 2800, B0*7.5);  // Milky Way / Earth end
  disc(yTop + B0*5, 0xffaacc, 2500, B0*6.5);  // Andromeda end
}

function updateLabels() {
  camera.updateMatrixWorld();
  const lMax = L_MAX * B0;
  function project(wx, wy, wz) {
    const v = new THREE.Vector3(wx, wy, wz).project(camera);
    return [(v.x+1)*.5*window.innerWidth, (1-(v.y+1)*.5)*window.innerHeight, v.z];
  }
  const [ex,ey,ez] = project(0, zL(-lMax,B0) - B0*5, 0);
  const [ax,ay,az] = project(0, zL( lMax,B0) + B0*5, 0);

  if (galaxiesOn && ez < 1) {
    lblEarth.style.display = 'block';
    lblEarth.style.left = (ex-55) + 'px';
    lblEarth.style.top  = (ey-18) + 'px';
  } else {
    lblEarth.style.display = 'none';
  }
  if (galaxiesOn && az < 1) {
    lblAndro.style.display = 'block';
    lblAndro.style.left = (ax-78) + 'px';
    lblAndro.style.top  = (ay-6) + 'px';
  } else {
    lblAndro.style.display = 'none';
  }
}

// ════════════════════════════════════════════════════════════════════
//  LIVE STATS
// ════════════════════════════════════════════════════════════════════
function updateStats() {
  // ρ(l=0) = −1/(8π b₀²)  [geometrized units, G=c=1]
  const rhoGeom = -1 / (8 * Math.PI * B0 * B0);
  document.getElementById('svRho').textContent =
    `−1/(8πb₀²) ≈ ${rhoGeom.toExponential(2)} [G=c=1]`;

  // Tidal-safe speed: γ²v² c² b₀/b₀³ ξ ≤ g_lim
  // → v_max ≈ √(g·b₀_SI / (c²·ξ))  (assuming 1 scene unit = 10 m, ξ=2 m)
  const b0_m = B0 * 10, C = 3e8;
  const vmax_frac = Math.sqrt(9.81 * b0_m / (C*C * 2.0));
  document.getElementById('svVmax').textContent =
    `${(vmax_frac*C).toFixed(2)} m/s  (1 unit=10 m, ξ=2 m)`;
}

// ════════════════════════════════════════════════════════════════════
//  REBUILD ALL on parameter change
// ════════════════════════════════════════════════════════════════════
function rebuild() {
  buildWormhole(); buildThroatRing(); buildRays();
  buildTraveler(); buildGalaxies(); updateStats();
}
rebuild();  // initial build

// ════════════════════════════════════════════════════════════════════
//  UI HANDLERS
// ════════════════════════════════════════════════════════════════════
document.getElementById('sB0').addEventListener('input', function() {
  B0 = parseFloat(this.value);
  document.getElementById('vB0').textContent = B0.toFixed(1);
  rebuild();
});
document.getElementById('sLen').addEventListener('input', function() {
  L_MAX = parseFloat(this.value);
  document.getElementById('vLen').textContent = L_MAX.toFixed(1) + '×b₀';
  rebuild();
});
document.getElementById('sSpd').addEventListener('input', function() {
  document.getElementById('vSpd').textContent = parseFloat(this.value).toFixed(1) + '×';
});

const btnWire   = document.getElementById('btnWire');
const btnRays   = document.getElementById('btnRays');
const btnTravel = document.getElementById('btnTravel');
const btnGalaxy = document.getElementById('btnGalaxy');

btnWire.addEventListener('click', () => {
  wireMode = !wireMode;
  btnWire.textContent = `Wireframe ${wireMode?'ON':'OFF'}`;
  btnWire.classList.toggle('act', wireMode);
  if (wormMesh) wormMesh.material.opacity = wireMode ? 0.10 : 0.66;
  if (wireMesh) wireMesh.visible = wireMode;
});

btnRays.addEventListener('click', () => {
  raysOn = !raysOn;
  btnRays.textContent = `Light Rays ${raysOn?'ON':'OFF'}`;
  btnRays.classList.toggle('off', !raysOn);
  rayLines.forEach(l => l.visible = raysOn);
});

btnTravel.addEventListener('click', () => {
  travelerOn = !travelerOn;
  btnTravel.textContent = `Traveler ${travelerOn?'ON':'OFF'}`;
  btnTravel.classList.toggle('off', !travelerOn);
  if (travelerMesh)  travelerMesh.visible  = travelerOn;
  if (travelerTrail) travelerTrail.visible = travelerOn;
});

btnGalaxy.addEventListener('click', () => {
  galaxiesOn = !galaxiesOn;
  btnGalaxy.textContent = `Galaxies ${galaxiesOn?'ON':'OFF'}`;
  btnGalaxy.classList.toggle('off', !galaxiesOn);
  galaxyObjs.forEach(g => g.visible = galaxiesOn);
});

// ════════════════════════════════════════════════════════════════════
//  ORBIT CONTROLS  (no external library needed)
// ════════════════════════════════════════════════════════════════════
let isDrag = false, isRight = false, lmx = 0, lmy = 0;

canvas.addEventListener('mousedown', e => {
  isDrag = true; isRight = e.button===2;
  lmx = e.clientX; lmy = e.clientY; e.preventDefault();
});
canvas.addEventListener('contextmenu', e => e.preventDefault());
window.addEventListener('mousemove', e => {
  if (!isDrag) return;
  const dx = e.clientX-lmx, dy = e.clientY-lmy;
  lmx = e.clientX; lmy = e.clientY;
  if (isRight) { panX -= dx*0.014; panY += dy*0.014; }
  else {
    camTheta += dx*0.007;
    camPhi = Math.max(0.05, Math.min(Math.PI-0.05, camPhi - dy*0.007));
  }
  updateCamera();
});
window.addEventListener('mouseup', () => isDrag = false);
canvas.addEventListener('wheel', e => {
  camRadius = Math.max(2, Math.min(200, camRadius * (1 + e.deltaY*0.0008)));
  updateCamera();
}, { passive:true });

// Touch support
let ltxy = null, lpd = null;
canvas.addEventListener('touchstart', e => {
  if (e.touches.length===1) ltxy={x:e.touches[0].clientX, y:e.touches[0].clientY};
  else if (e.touches.length===2)
    lpd = Math.hypot(e.touches[0].clientX-e.touches[1].clientX,
                     e.touches[0].clientY-e.touches[1].clientY);
  e.preventDefault();
}, { passive:false });
canvas.addEventListener('touchmove', e => {
  if (e.touches.length===1 && ltxy) {
    const dx=e.touches[0].clientX-ltxy.x, dy=e.touches[0].clientY-ltxy.y;
    ltxy={x:e.touches[0].clientX, y:e.touches[0].clientY};
    camTheta += dx*0.007;
    camPhi = Math.max(0.05, Math.min(Math.PI-0.05, camPhi - dy*0.007));
    updateCamera();
  } else if (e.touches.length===2 && lpd) {
    const d=Math.hypot(e.touches[0].clientX-e.touches[1].clientX,
                       e.touches[0].clientY-e.touches[1].clientY);
    camRadius = Math.max(2, Math.min(200, camRadius*lpd/d));
    lpd = d; updateCamera();
  }
  e.preventDefault();
}, { passive:false });
canvas.addEventListener('touchend', () => { ltxy=null; lpd=null; });

window.addEventListener('resize', () => {
  renderer.setSize(window.innerWidth, window.innerHeight);
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
});

// ════════════════════════════════════════════════════════════════════
//  ANIMATION LOOP
// ════════════════════════════════════════════════════════════════════
let lastT = performance.now();

function animate() {
  requestAnimationFrame(animate);
  const now = performance.now();
  const dt  = Math.min((now - lastT) * 0.001, 0.05);
  lastT = now;
  const speed = parseFloat(document.getElementById('sSpd').value);

  // Pulse throat glow
  throatGlow.intensity = 2.8 + 2.0*Math.sin(now*0.0024);

  // Slowly spin galaxy discs
  for (const g of galaxyObjs) g.rotation.y += dt * 0.055;

  // Animate traveler
  if (travelerOn && travelerMesh) stepTraveler(dt, speed);

  updateLabels();
  renderer.render(scene, camera);
}
animate();
</script>
</body>
</html>
