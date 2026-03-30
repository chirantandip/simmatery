// ============================================================================
// SimMatery — Physics Simulation Engine
// ============================================================================
//
// ARCHITECTURE OVERVIEW
//
//   This file uses a model-config-dispatch pattern. Each physics model is
//   defined by 3 functions (init, step, draw) registered in modelConfigs.
//   Global function references (sim_initFn, sim_stepFn, sim_drawFn) are set
//   when a model is selected, eliminating all if/else branching from the
//   simulation loop.
//
//   The file is divided into two types of sections:
//
//   [DO NOT MODIFY]  Framework infrastructure — the simulation engine,
//                    rendering pipeline, UI system, and event wiring.
//                    Changing these can break ALL models.
//
//   [EXTENSIBLE]     Model-specific code — add your new physics model here.
//                    These sections are safe to extend without affecting
//                    existing models.
//
// ============================================================================
// HOW TO ADD A NEW PHYSICS MODEL — 6 STEP CHECKLIST
// ============================================================================
//
//   Step 1 — Declare parameter variables (Section: Parameter Variables)
//            Add global variables for your model's tuneable parameters.
//
//   Step 2 — Write an init function (Section: Model Init Functions)
//            Allocate your fields (Float32Array) and set initial conditions.
//            Signature: () => void. Access GridSizeX, GridSizeY from globals.
//
//   Step 3 — Write a step function (Section: Model Step Functions)
//            Implement one simulation tick. Signature: () => void.
//            Read/write your fields in-place or swap with _NEXT buffers.
//
//   Step 4 — Write or choose a draw function (Section: Model Draw Functions)
//            Map your field data to pixel colors via the colorPalette.
//            If your field is ∈ [0,1], reuse drawField01. If ±1, reuse
//            drawFieldIsing. Otherwise write a new one.
//
//   Step 5 — Register in modelConfigs (Section: Model Config)
//            Add one line: yourModelKey: { init, step, draw }
//            The key MUST match data-model in index.html dropdown.
//
//   Step 6 — Add paramConfigs entry (Section: Param Configs)
//            Add a slider config array for your model's UI parameters.
//
//   Also add:
//     - A dropdown item in index.html: <a data-model="yourkey">...</a>
//     - Your model's display name to displayWelcomeCanvas models list
//
// ============================================================================


// ============================================================================
// CANVAS SETUP — DO NOT MODIFY
// ============================================================================

let simCTX  = document.getElementById('canvas').getContext('2d', { alpha: false });
let canvas  = simCTX.canvas;

let simModel = 'cahnhilliard';
let simRequest;
let simRunning = false;
let lastTime = 0;

let cellSize = 8;

// ============================================================================
// MODEL KEYS — DO NOT MODIFY (add yours here, keep existing ones)
// ============================================================================

const simCahn  = 'cahnhilliard';
const simIsing = 'ising';
const simDLA   = 'dla';
const simGray  = 'gray';
const simAllen = 'allencahn';
const simKob   = 'kobayashi';

// ============================================================================
// FIELD DECLARATIONS — EXTENSIBLE
// ============================================================================
// Add your model-specific field arrays here. Declare them at this level so
// they persist across frames. Use Float32Array for grid data.
//
//   FIELD       — shared by Cahn-Hilliard, Allen-Cahn, Ising, DLA
//   FIELD_NEXT  — swap buffer for FIELD (Cahn-Hilliard, Allen-Cahn)
//   FIELD_MU    — chemical potential (Cahn-Hilliard only)
//   U_FIELD, V_FIELD, U_NEXT, V_NEXT — dual-field (Gray-Scott only)
//
// If your model needs its own field arrays, declare them here.
// ============================================================================

/// generic
let FIELD;
let FIELD_NEXT;

/// for cahn hillard
let FIELD_MU;

/// for gray scott
let U_FIELD, U_NEXT;
let V_FIELD, V_NEXT;

// ============================================================================
// GRID DIMENSIONS — DO NOT MODIFY
// ============================================================================
// These are set by sim_resizeCanvas() based on container size and cellSize.
// All model code should read GridSizeX / GridSizeY from these globals.
// ============================================================================

let GridSizeX = 128;
let GridSizeY = 128;

// ============================================================================
// PARAMETER VARIABLES — EXTENSIBLE
// ============================================================================
// Add your model's tuneable parameters here as global `let` variables.
// These are read by your step function and bound to UI sliders via
// paramConfigs (see Section: Param Configs).
//
// Naming convention: modelPrefix_paramName (e.g., cahn_M, ising_T, ac_F)
// ============================================================================

let cahn_dx = 1.0;
let cahn_dt = 0.25;

let cahn_M  = 0.25; /// Mobility
let cahn_K  = 0.25; /// Kappa
let cahn_W  = 0.99; /// Double well energy

let ising_T = 2.0;  /// Temperature
let ising_B = 0.0;  /// Magnetic field
let ising_J = 1.2;  /// Coupling constant
let ising_N = 1;    /// neighborhood

let dla_P   = 0.5;  /// Sticking probability
let dla_N   = 10;   /// Number of particles
let dla_D   = 5;    /// Launch distance from rightmost occupied cell

let dla_Rightmost = 0;

let gs_Du   = 0.16;
let gs_Dv   = 0.08;
let gs_F    = 0.025;
let gs_k    = 0.06;

let ac_M    = 1.0;   /// Mobility
let ac_K    = 0.5;   /// Interfacial energy (kappa)
let ac_W    = 0.9;  /// Double-well coefficient
let ac_F    = 0.0;   /// Driving force

let kobayashi_epsilon = 1.0;   // Interface thickness parameter
let kobayashi_mobility = 1.0;  // Mobility
let kobayashi_undercooling = -0.1;  // Constant undercooling
let kobayashi_dt = 0.3;        // Time step

// ============================================================================
// CONSTANTS — DO NOT MODIFY
// ============================================================================

const barWidth = 400;
const cellDim  = 3;
const cartDirections = [ [0, 1], [0, -1], [1, 0], [-1, 0] ];

// ============================================================================
// FIELD ALLOCATION — DO NOT MODIFY
// ============================================================================
// Conditional allocation: reuses existing Float32Array if size matches,
// allocates new only when size changes or array is null. Preserves field
// values across model switches (where compatible).
// ============================================================================

let _prevModel = null;

function allocIfNeeded(arr, size) {
    if (!arr || arr.length !== size) {
        return new Float32Array(size);
    }
    return arr;
}

// ============================================================================
// GRID LIFECYCLE — DO NOT MODIFY
// ============================================================================
// sim_resizeCanvas: resizes canvas to fill container, recomputes grid size
// sim_initModelGrid: dispatches to the current model's init function
// ============================================================================

function sim_GetSizeY(SizeX) {
    return Math.floor(canvas.height / (canvas.width / SizeX) + 1);
}

function sim_resizeCanvas() {
    const container = canvas.parentElement;
    let w = container.clientWidth;
    let h = container.clientHeight;

    // Ensure minimum size
    w = Math.max(w, 100);
    h = Math.max(h, 100);

    // Canvas fills entire container
    canvas.width = w;
    canvas.height = h;

    // Simulation grid = canvas / cellSize
    GridSizeX = Math.floor(w / cellSize);
    GridSizeY = Math.floor(h / cellSize);

    // Ensure minimum simulation grid
    GridSizeX = Math.max(GridSizeX, 16);
    GridSizeY = Math.max(GridSizeY, 16);

    // Reset ImageData to force recreation
    simImageData = null;

    sim_initModelGrid();
    sim_drawGrid();
}

// CRITICAL: sim_initModelGrid is the init dispatcher. Do not add if/else here.
// The current model's init function is called via sim_initFn.
function sim_initModelGrid() {
    sim_Stop();
    sim_initFn();
}

// ============================================================================
// COLOR PALETTE — DO NOT MODIFY
// ============================================================================
// Linear interpolation between two RGB colors. All draw functions use this.
// ============================================================================

let colorPalette = new Uint8Array(256 * 3);
let simImageData = null;

function initColorPalette() {
    const rgbA = [13, 94, 166];
    const rgbB = [234, 166, 77];
    for (let i = 0; i < 256; i++) {
        const t = i / 255;
        colorPalette[i * 3 + 0] = Math.floor(rgbA[0] + (rgbB[0] - rgbA[0]) * t);
        colorPalette[i * 3 + 1] = Math.floor(rgbA[1] + (rgbB[1] - rgbA[1]) * t);
        colorPalette[i * 3 + 2] = Math.floor(rgbA[2] + (rgbB[2] - rgbA[2]) * t);
    }
}

// ============================================================================
// WELCOME SCREEN — DO NOT MODIFY (except: add your model name to the list)
// ============================================================================

const USE_MAXIMALIST_WELCOME = true;

let welcomeAnimFrame = 0;
let stars = [];
let mangoes = [];
let crystals = [];
let shootingStars = [];
let atoms = [];
let butterflies = [];
let sparkles = [];
let comets = [];
let welcomeRunning = false;
let welcomeAnimId = null;

function initWelcomeParticles() {
    stars = [];
    mangoes = [];
    crystals = [];
    shootingStars = [];
    atoms = [];
    butterflies = [];
    sparkles = [];
    comets = [];
    
    for (let i = 0; i < 80; i++) {
        stars.push({
            x: Math.random() * canvas.width,
            y: Math.random() * canvas.height,
            size: 1 + Math.random() * 2,
            color: ['#FFFFFF', '#00FFFF', '#FFFF00', '#FF69B4'][Math.floor(Math.random() * 4)],
            phase: Math.random() * Math.PI * 2,
            speed: 0.02 + Math.random() * 0.03
        });
    }
    
    for (let i = 0; i < 12; i++) {
        const side = Math.floor(Math.random() * 4);
        let x, y, speedX, speedY;
        
        if (side === 0) { // right
            x = canvas.width + 50;
            y = Math.random() * canvas.height;
            speedX = -(0.3 + Math.random() * 0.5);
            speedY = (Math.random() - 0.5) * 0.3;
        } else if (side === 1) { // left
            x = -50;
            y = Math.random() * canvas.height;
            speedX = 0.3 + Math.random() * 0.5;
            speedY = (Math.random() - 0.5) * 0.3;
        } else if (side === 2) { // top
            x = Math.random() * canvas.width;
            y = -50;
            speedX = (Math.random() - 0.5) * 0.3;
            speedY = 0.3 + Math.random() * 0.5;
        } else { // bottom
            x = Math.random() * canvas.width;
            y = canvas.height + 50;
            speedX = (Math.random() - 0.5) * 0.3;
            speedY = -(0.3 + Math.random() * 0.5);
        }
        
        mangoes.push({
            x: x,
            y: y,
            size: 18 + Math.random() * 12,
            speedX: speedX,
            speedY: speedY,
            phase: Math.random() * Math.PI * 2,
            amplitude: 20 + Math.random() * 40
        });
    }
    
    for (let i = 0; i < 10; i++) {
        const side = Math.floor(Math.random() * 4);
        let x, y, speedX, speedY;
        
        if (side === 0) { // right
            x = canvas.width + 50;
            y = Math.random() * canvas.height;
            speedX = -(0.3 + Math.random() * 0.4);
            speedY = (Math.random() - 0.5) * 0.2;
        } else if (side === 1) { // left
            x = -50;
            y = Math.random() * canvas.height;
            speedX = 0.3 + Math.random() * 0.4;
            speedY = (Math.random() - 0.5) * 0.2;
        } else if (side === 2) { // top
            x = Math.random() * canvas.width;
            y = -50;
            speedX = (Math.random() - 0.5) * 0.2;
            speedY = 0.4 + Math.random() * 0.4;
        } else { // bottom
            x = Math.random() * canvas.width;
            y = canvas.height + 50;
            speedX = (Math.random() - 0.5) * 0.2;
            speedY = -(0.4 + Math.random() * 0.4);
        }
        
        crystals.push({
            x: x,
            y: y,
            size: 14 + Math.random() * 10,
            rotation: Math.random() * Math.PI * 2,
            speedX: speedX,
            speedY: speedY,
            rotationSpeed: 0.01 + Math.random() * 0.02,
            color: ['#00FFFF', '#FF00FF', '#FFD700', '#39FF14'][Math.floor(Math.random() * 4)]
        });
    }
    
    for (let i = 0; i < 4; i++) {
        shootingStars.push({
            x: Math.random() * canvas.width,
            y: Math.random() * canvas.height * 0.6,
            length: 30 + Math.random() * 40,
            speed: 4 + Math.random() * 3,
            angle: Math.PI / 4 + (Math.random() - 0.5) * 0.3,
            active: Math.random() > 0.7,
            cooldown: Math.random() * 200
        });
    }
    
    for (let i = 0; i < 20; i++) {
        atoms.push({
            x: Math.random() * canvas.width,
            y: Math.random() * canvas.height,
            size: 10 + Math.random() * 8,
            rotation: Math.random() * Math.PI * 2,
            rotationSpeed: 0.02 + Math.random() * 0.03,
            orbitRadius: 12 + Math.random() * 8,
            orbitTilt: Math.random() * Math.PI,
            color: ['#00FFFF', '#FF00FF', '#FFD700', '#39FF14'][Math.floor(Math.random() * 4)]
        });
    }
    
    for (let i = 0; i < 4; i++) {
        butterflies.push({
            x: Math.random() > 0.5 ? -30 : canvas.width + 30,
            y: Math.random() * canvas.height,
            size: 12 + Math.random() * 8,
            speedX: (Math.random() > 0.5 ? 1 : -1) * (0.3 + Math.random() * 0.4),
            phase: Math.random() * Math.PI * 2,
            flutterSpeed: 0.15 + Math.random() * 0.1,
            color: ['#FF1493', '#FF6600', '#FFD700', '#FF69B4'][Math.floor(Math.random() * 4)]
        });
    }
    
    for (let i = 0; i < 6; i++) {
        sparkles.push({
            x: Math.random() * canvas.width,
            y: Math.random() * canvas.height,
            size: 4 + Math.random() * 6,
            life: 0,
            maxLife: 30 + Math.random() * 40,
            delay: Math.random() * 150,
            color: ['#FFFFFF', '#FFFF00', '#00FFFF', '#FF69B4'][Math.floor(Math.random() * 4)]
        });
    }
    
    for (let i = 0; i < 12; i++) {
        comets.push({
            x: Math.random() * canvas.width,
            y: Math.random() * canvas.height * 0.5,
            size: 6 + Math.random() * 4,
            speed: 1.5 + Math.random() * 1,
            angle: Math.PI / 4 + (Math.random() - 0.5) * 0.5,
            trailLength: 60 + Math.random() * 40,
            active: Math.random() > 0.6,
            cooldown: Math.random() * 300,
            color: ['#00FFFF', '#FF00FF', '#FFD700', '#FF69B4', '#39FF14'][Math.floor(Math.random() * 5)]
        });
    }
}

function drawStar(s) {
    const twinkle = 0.4 + 0.6 * Math.abs(Math.sin(welcomeAnimFrame * s.speed + s.phase));
    simCTX.save();
    simCTX.globalAlpha = twinkle;
    simCTX.fillStyle = s.color;
    simCTX.beginPath();
    simCTX.arc(s.x, s.y, s.size, 0, Math.PI * 2);
    simCTX.fill();
    simCTX.restore();
}

function drawMango(m) {
    simCTX.save();
    simCTX.translate(m.x, m.y);
    
    const gradient = simCTX.createLinearGradient(-m.size/2, 0, m.size/2, 0);
    gradient.addColorStop(0, '#FFD700');
    gradient.addColorStop(0.5, '#FF8C00');
    gradient.addColorStop(1, '#FFD700');
    simCTX.fillStyle = gradient;
    
    simCTX.beginPath();
    simCTX.ellipse(0, 0, m.size * 0.6, m.size, 0, 0, Math.PI * 2);
    simCTX.fill();
    
    simCTX.fillStyle = '#228B22';
    simCTX.beginPath();
    simCTX.moveTo(0, -m.size);
    simCTX.lineTo(-3, -m.size - 8);
    simCTX.lineTo(3, -m.size - 8);
    simCTX.closePath();
    simCTX.fill();
    
    simCTX.restore();
}

function drawCrystal(c) {
    simCTX.save();
    simCTX.translate(c.x, c.y);
    simCTX.rotate(c.rotation);
    
    simCTX.strokeStyle = c.color;
    simCTX.lineWidth = 2;
    simCTX.shadowColor = c.color;
    simCTX.shadowBlur = 10;
    
    const s = c.size;
    simCTX.beginPath();
    simCTX.moveTo(0, -s);
    simCTX.lineTo(s * 0.5, 0);
    simCTX.lineTo(0, s);
    simCTX.lineTo(-s * 0.5, 0);
    simCTX.closePath();
    simCTX.stroke();
    
    simCTX.fillStyle = c.color + '44';
    simCTX.fill();
    
    simCTX.restore();
}

function drawShootingStar(s) {
    if (!s.active) return;
    
    const endX = s.x - Math.cos(s.angle) * s.length;
    const endY = s.y - Math.sin(s.angle) * s.length;
    
    const gradient = simCTX.createLinearGradient(s.x, s.y, endX, endY);
    gradient.addColorStop(0, '#FFFFFF');
    gradient.addColorStop(0.3, '#FFFF00');
    gradient.addColorStop(1, 'transparent');
    
    simCTX.save();
    simCTX.strokeStyle = gradient;
    simCTX.lineWidth = 2;
    simCTX.lineCap = 'round';
    simCTX.shadowColor = '#FFFF00';
    simCTX.shadowBlur = 10;
    
    simCTX.beginPath();
    simCTX.moveTo(s.x, s.y);
    simCTX.lineTo(endX, endY);
    simCTX.stroke();
    
    simCTX.fillStyle = '#FFFFFF';
    simCTX.beginPath();
    simCTX.arc(s.x, s.y, 2, 0, Math.PI * 2);
    simCTX.fill();
    
    simCTX.restore();
}

function drawAtom(a) {
    simCTX.save();
    simCTX.translate(a.x, a.y);
    
    simCTX.fillStyle = a.color;
    simCTX.shadowColor = a.color;
    simCTX.shadowBlur = 8;
    simCTX.beginPath();
    simCTX.arc(0, 0, a.size * 0.4, 0, Math.PI * 2);
    simCTX.fill();
    
    const orbitX = Math.cos(a.rotation) * a.orbitRadius;
    const orbitY = Math.sin(a.rotation) * a.orbitRadius * 0.5;
    
    simCTX.strokeStyle = a.color + '88';
    simCTX.lineWidth = 1.5;
    simCTX.beginPath();
    simCTX.ellipse(0, 0, a.orbitRadius, a.orbitRadius * 0.5, a.orbitTilt, 0, Math.PI * 2);
    simCTX.stroke();
    
    simCTX.fillStyle = '#FFFFFF';
    simCTX.beginPath();
    simCTX.arc(orbitX, orbitY, 3, 0, Math.PI * 2);
    simCTX.fill();
    
    const orbit2X = Math.cos(a.rotation + Math.PI) * a.orbitRadius;
    const orbit2Y = Math.sin(a.rotation + Math.PI) * a.orbitRadius * 0.5;
    
    simCTX.beginPath();
    simCTX.arc(orbit2X, orbit2Y, 3, 0, Math.PI * 2);
    simCTX.fill();
    
    simCTX.restore();
}

function drawButterfly(b) {
    simCTX.save();
    simCTX.translate(b.x, b.y);
    
    const flutter = Math.sin(welcomeAnimFrame * b.flutterSpeed * 10 + b.phase) * 0.4;
    simCTX.rotate(flutter);
    
    simCTX.fillStyle = b.color;
    simCTX.shadowColor = b.color;
    simCTX.shadowBlur = 8;
    
    simCTX.beginPath();
    simCTX.ellipse(-b.size * 0.4, 0, b.size * 0.5, b.size * 0.3, 0, 0, Math.PI * 2);
    simCTX.fill();
    simCTX.beginPath();
    simCTX.ellipse(b.size * 0.4, 0, b.size * 0.5, b.size * 0.3, 0, 0, Math.PI * 2);
    simCTX.fill();
    
    simCTX.strokeStyle = '#000000';
    simCTX.lineWidth = 1;
    simCTX.beginPath();
    simCTX.moveTo(0, -b.size * 0.4);
    simCTX.lineTo(0, b.size * 0.4);
    simCTX.stroke();
    
    simCTX.restore();
}

function drawSparkle(s) {
    if (s.delay > 0) return;
    
    const lifeRatio = s.life / s.maxLife;
    if (lifeRatio > 1) return;
    
    const alpha = lifeRatio < 0.2 ? lifeRatio * 5 : (1 - lifeRatio);
    const size = s.size * (lifeRatio < 0.2 ? lifeRatio * 5 : (1 - lifeRatio * 0.5));
    
    simCTX.save();
    simCTX.globalAlpha = alpha;
    simCTX.fillStyle = s.color;
    simCTX.shadowColor = s.color;
    simCTX.shadowBlur = 10;
    
    const half = size / 2;
    simCTX.beginPath();
    simCTX.moveTo(0, -size);
    simCTX.lineTo(half * 0.3, -half * 0.3);
    simCTX.lineTo(size, 0);
    simCTX.lineTo(half * 0.3, half * 0.3);
    simCTX.lineTo(0, size);
    simCTX.lineTo(-half * 0.3, half * 0.3);
    simCTX.lineTo(-size, 0);
    simCTX.lineTo(-half * 0.3, -half * 0.3);
    simCTX.closePath();
    simCTX.fill();
    
    simCTX.restore();
}

function drawComet(c) {
    if (!c.active) return;
    
    const trailGradient = simCTX.createLinearGradient(
        c.x, c.y,
        c.x - Math.cos(c.angle) * c.trailLength,
        c.y - Math.sin(c.angle) * c.trailLength
    );
    trailGradient.addColorStop(0, '#FFFFFF');
    trailGradient.addColorStop(0.2, c.color || '#00FFFF');
    trailGradient.addColorStop(1, 'transparent');
    
    simCTX.save();
    simCTX.strokeStyle = trailGradient;
    simCTX.lineWidth = c.size;
    simCTX.lineCap = 'round';
    simCTX.shadowColor = '#FFFFFF';
    simCTX.shadowBlur = 15;
    
    simCTX.beginPath();
    simCTX.moveTo(c.x, c.y);
    simCTX.lineTo(
        c.x - Math.cos(c.angle) * c.trailLength,
        c.y - Math.sin(c.angle) * c.trailLength
    );
    simCTX.stroke();
    
    simCTX.fillStyle = '#FFFFFF';
    simCTX.beginPath();
    simCTX.arc(c.x, c.y, c.size * 0.6, 0, Math.PI * 2);
    simCTX.fill();
    
    simCTX.restore();
}

function drawMaximalistWelcome() {
    welcomeAnimFrame += 0.016;
    
    const gradient = simCTX.createLinearGradient(0, 0, canvas.width, canvas.height);
    gradient.addColorStop(0, '#8B0000');
    gradient.addColorStop(0.4, '#FF1493');
    gradient.addColorStop(0.7, '#FF6600');
    gradient.addColorStop(1, '#FFD700');
    simCTX.fillStyle = gradient;
    simCTX.fillRect(0, 0, canvas.width, canvas.height);
    
    stars.forEach(drawStar);
    
    mangoes.forEach(m => {
        m.x += m.speedX;
        m.y += m.speedY + Math.sin(welcomeAnimFrame * 2 + m.phase) * m.amplitude * 0.02;
        
        const margin = 60;
        if (m.x < -margin || m.x > canvas.width + margin || 
            m.y < -margin || m.y > canvas.height + margin) {
            const side = Math.floor(Math.random() * 4);
            if (side === 0) { m.x = canvas.width + 30; m.y = Math.random() * canvas.height; }
            else if (side === 1) { m.x = -30; m.y = Math.random() * canvas.height; }
            else if (side === 2) { m.x = Math.random() * canvas.width; m.y = -30; }
            else { m.x = Math.random() * canvas.width; m.y = canvas.height + 30; }
        }
        drawMango(m);
    });
    
    crystals.forEach(c => {
        c.x += c.speedX;
        c.y += c.speedY;
        c.rotation += c.rotationSpeed;
        
        const margin = 60;
        if (c.x < -margin || c.x > canvas.width + margin || 
            c.y < -margin || c.y > canvas.height + margin) {
            const side = Math.floor(Math.random() * 4);
            if (side === 0) { c.x = canvas.width + 30; c.y = Math.random() * canvas.height; }
            else if (side === 1) { c.x = -30; c.y = Math.random() * canvas.height; }
            else if (side === 2) { c.x = Math.random() * canvas.width; c.y = -30; }
            else { c.x = Math.random() * canvas.width; c.y = canvas.height + 30; }
        }
        drawCrystal(c);
    });
    
    shootingStars.forEach(s => {
        if (!s.active) {
            s.cooldown--;
            if (s.cooldown <= 0) {
                s.active = true;
                s.x = canvas.width + 50;
                s.y = Math.random() * canvas.height * 0.5;
            }
        } else {
            s.x -= Math.cos(s.angle) * s.speed;
            s.y -= Math.sin(s.angle) * s.speed;
            if (s.x < -100 || s.y < -100) {
                s.active = false;
                s.cooldown = 150 + Math.random() * 200;
            }
        }
        drawShootingStar(s);
    });
    
    atoms.forEach(a => {
        a.rotation += a.rotationSpeed;
        drawAtom(a);
    });
    
    butterflies.forEach(b => {
        b.x += b.speedX;
        b.y += Math.sin(welcomeAnimFrame * 2 + b.phase) * 0.5;
        if ((b.speedX > 0 && b.x > canvas.width + 50) || (b.speedX < 0 && b.x < -50)) {
            b.x = b.speedX > 0 ? -30 : canvas.width + 30;
            b.y = Math.random() * canvas.height;
        }
        drawButterfly(b);
    });
    
    sparkles.forEach(s => {
        if (s.delay > 0) {
            s.delay--;
        } else {
            s.life++;
            if (s.life > s.maxLife) {
                s.x = Math.random() * canvas.width;
                s.y = Math.random() * canvas.height;
                s.life = 0;
                s.maxLife = 30 + Math.random() * 40;
            }
        }
        drawSparkle(s);
    });
    
    comets.forEach(c => {
        if (!c.active) {
            c.cooldown--;
            if (c.cooldown <= 0) {
                c.active = true;
                c.x = -50;
                c.y = Math.random() * canvas.height * 0.4;
                c.color = ['#00FFFF', '#FF00FF', '#FFD700', '#FF69B4', '#39FF14'][Math.floor(Math.random() * 5)];
            }
        } else {
            c.x += Math.cos(c.angle) * c.speed;
            c.y += Math.sin(c.angle) * c.speed;
            if (c.x > canvas.width + 100 || c.y > canvas.height + 100) {
                c.active = false;
                c.cooldown = 200 + Math.random() * 300;
            }
        }
        drawComet(c);
    });
    
    const cy = canvas.height / 2;
    const cx = canvas.width / 2;
    
    const glowCycle = (Math.floor(welcomeAnimFrame * 2) % 4);
    const glowColors = ['#00FFFF', '#FF00FF', '#FFD700', '#39FF14'];
    const currentGlow = glowColors[glowCycle];
    
    const bobOffset = Math.sin(welcomeAnimFrame * 3) * 3;
    
    simCTX.save();
    simCTX.shadowColor = '#000000';
    simCTX.shadowBlur = 0;
    simCTX.shadowOffsetX = 4;
    simCTX.shadowOffsetY = 4;
    simCTX.fillStyle = '#000000';
    simCTX.font = '100px "ComputerModernSans", sans-serif';
    simCTX.textAlign = 'center';
    simCTX.fillText('SimMatery', cx + 4, cy - 50 + bobOffset + 4);
    simCTX.restore();
    
    simCTX.save();
    simCTX.shadowColor = '#FF00FF';
    simCTX.shadowBlur = 15;
    simCTX.fillStyle = '#FF00FF';
    simCTX.font = '100px "ComputerModernSans", sans-serif';
    simCTX.textAlign = 'center';
    simCTX.fillText('SimMatery', cx - 3, cy - 50 + bobOffset - 3);
    simCTX.restore();
    
    simCTX.save();
    simCTX.shadowColor = '#00FFFF';
    simCTX.shadowBlur = 15;
    simCTX.fillStyle = '#00FFFF';
    simCTX.font = '100px "ComputerModernSans", sans-serif';
    simCTX.textAlign = 'center';
    simCTX.fillText('SimMatery', cx + 3, cy - 50 + bobOffset + 3);
    simCTX.restore();
    
    simCTX.save();
    simCTX.shadowColor = currentGlow;
    simCTX.shadowBlur = 25;
    const titleGradient = simCTX.createLinearGradient(cx - 150, 0, cx + 150, 0);
    titleGradient.addColorStop(0, '#FF1493');
    titleGradient.addColorStop(0.5, '#FFD700');
    titleGradient.addColorStop(1, '#00FFFF');
    simCTX.fillStyle = titleGradient;
    simCTX.font = '100px "ComputerModernSans", sans-serif';
    simCTX.textAlign = 'center';
    simCTX.fillText('SimMatery', cx, cy - 50 + bobOffset);
    simCTX.restore();
    
    const time = welcomeAnimFrame;
    simCTX.strokeStyle = '#FFD700';
    simCTX.lineWidth = 3;
    simCTX.shadowColor = '#FFD700';
    simCTX.shadowBlur = 10;
    
    const borderY = 30;
    simCTX.beginPath();
    for (let x = 0; x <= canvas.width; x += 5) {
        const y = borderY + Math.sin(x * 0.02 + time * 2) * 8;
        if (x === 0) simCTX.moveTo(x, y);
        else simCTX.lineTo(x, y);
    }
    simCTX.stroke();
    
    const bottomBorderY = canvas.height - 30;
    simCTX.beginPath();
    for (let x = 0; x <= canvas.width; x += 5) {
        const y = bottomBorderY + Math.sin(x * 0.02 + time * 2 + Math.PI) * 8;
        if (x === 0) simCTX.moveTo(x, y);
        else simCTX.lineTo(x, y);
    }
    simCTX.stroke();
    
    if (!welcomeRunning) {
        welcomeRunning = true;
        animateWelcome();
    }
}

function animateWelcome() {
    if (!welcomeRunning) {
        welcomeAnimId = null;
        return;
    }
    displayWelcomeCanvas();
    welcomeAnimId = requestAnimationFrame(animateWelcome);
}

function displayWelcomeCanvas() {
    if (USE_MAXIMALIST_WELCOME) {
        drawMaximalistWelcome();
    } else {
        drawSimpleWelcome();
    }
}

function drawSimpleWelcome() {
    simCTX.fillStyle = '#0f172a';
    simCTX.fillRect(0, 0, canvas.width, canvas.height);

    const cx = canvas.width / 2;
    const cy = canvas.height / 2;

    simCTX.fillStyle = '#60a5fa';
    simCTX.font = 'bold 48px -apple-system, BlinkMacSystemFont, sans-serif';
    simCTX.textAlign = 'center';
    simCTX.fillText('SimMatery', cx, cy - 50);

    simCTX.fillStyle = '#d1d5db';
    simCTX.font = '18px -apple-system, BlinkMacSystemFont, sans-serif';
    simCTX.fillText('Select a model from the dropdown to begin', cx, cy + 20);

    simCTX.fillStyle = '#6b7280';
    simCTX.font = '14px -apple-system, BlinkMacSystemFont, sans-serif';
     // EXTENSIBLE: add your model display name to this list
     const models = ['Cahn-Hilliard Model', 'Allen-Cahn Model', 'Ising Model', 'DLA', 'Gray-Scott', 'Kobayashi Model'];
     models.forEach((m, i) => simCTX.fillText('• ' + m, cx, cy + 60 + i * 20));
}

// ============================================================================
// DRAW DISPATCH — DO NOT MODIFY
// ============================================================================
// CRITICAL: sim_drawGrid dispatches to the current model's draw function.
// Do not add if/else here. The current draw function is called via sim_drawFn.
// ============================================================================

function sim_drawGrid() {
    sim_drawFn();
}

// ============================================================================
// SHARED GRID UTILITIES — DO NOT MODIFY
// ============================================================================
// These are used by step and draw functions across all models.
// ============================================================================

function idx(x, y) {
    return ((x + GridSizeX) % GridSizeX) + ((y + GridSizeY) % GridSizeY) * GridSizeX;
}

// 9-point stencil Laplacian with toroidal (wraparound) boundary conditions.
// Used by: Cahn-Hilliard, Allen-Cahn, Gray-Scott.
// Signature: lapGrid(fieldArray, x, y) => number
function lapGrid(field, x, y) {
    const w = GridSizeX;
    const left = field[((x - 1 + w) % w) + y * w];
    const right = field[((x + 1) % w) + y * w];
    const up = field[x + ((y - 1 + GridSizeY) % GridSizeY) * w];
    const down = field[x + ((y + 1 + GridSizeY) % GridSizeY) * w];
    const ul = field[((x - 1 + w) % w) + ((y - 1 + GridSizeY) % GridSizeY) * w];
    const ur = field[((x + 1) % w) + ((y - 1 + GridSizeY) % GridSizeY) * w];
    const dl = field[((x - 1 + w) % w) + ((y + 1 + GridSizeY) % GridSizeY) * w];
    const dr = field[((x + 1) % w) + ((y + 1 + GridSizeY) % GridSizeY) * w];
    const center = field[x + y * w];
    return (4 * (left + right + up + down) + (ul + ur + dl + dr) - 20 * center) / (6 * cahn_dx * cahn_dx);
}


// ============================================================================
// ============================================================================
// ============================================================================
// MODEL STEP FUNCTIONS — EXTENSIBLE
// ============================================================================
// ============================================================================
// ============================================================================
// Each model has a step function that runs one simulation tick.
//
// RULES:
//   - Signature: () => void  (no arguments, reads globals)
//   - Read from: FIELD, U_FIELD, V_FIELD (your model's input arrays)
//   - Write to:  FIELD_NEXT, U_NEXT, V_NEXT (swap buffers)
//   - Swap at end: [FIELD, FIELD_NEXT] = [FIELD_NEXT, FIELD]
//   - Use GridSizeX, GridSizeY for loop bounds
//   - Use lapGrid(field, x, y) for Laplacian if needed
//   - Read parameters from your global variables (e.g., cahn_M, ising_T)
//   - Do NOT call sim_drawGrid() or requestAnimationFrame — the run loop
//     handles that.
//
// To add a new step function, write it here and reference it in modelConfigs.
// ============================================================================

function sim_cahnHilliardStep() {
    const w = GridSizeX;
    const h = GridSizeY;

    // Compute chemical potential
    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            const i = x + y * w;
            const phi = FIELD[i];
            const dwdphi = 2 * cahn_W * phi * (1 - phi) * (1 - 2 * phi);
            FIELD_MU[i] = dwdphi - cahn_K * lapGrid(FIELD, x, y);
        }
    }
    // Update field
    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            const i = x + y * w;
            let newVal = FIELD[i] + cahn_dt * cahn_M * lapGrid(FIELD_MU, x, y);
            FIELD_NEXT[i] = newVal < 0 ? 0 : (newVal > 1 ? 1 : newVal);
        }
    }
    [FIELD, FIELD_NEXT] = [FIELD_NEXT, FIELD];
}

function sim_allenCahnStep() {
    const w = GridSizeX;
    const h = GridSizeY;

    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            const i = x + y * w;
            const phi = FIELD[i];
            const dwdphi = 2 * ac_W * phi * (1 - phi) * (1 - 2 * phi);
            const lap = lapGrid(FIELD, x, y);
            const dphi = -ac_M * (dwdphi - ac_K * lap - ac_F);
            let newVal = phi + cahn_dt * dphi;
            FIELD_NEXT[i] = newVal < 0 ? 0 : (newVal > 1 ? 1 : newVal);
        }
    }
    [FIELD, FIELD_NEXT] = [FIELD_NEXT, FIELD];
}

function sim_isingStep() {
    const w = GridSizeX;
    const h = GridSizeY;
    const total = w * h;

    for (let k = 0; k < total; k++) {
        const x = Math.floor(Math.random() * w);
        const y = Math.floor(Math.random() * h);
        const i = x + y * w;

        const spin = FIELD[i];
        let sum = 0;

        if (ising_N === 1) {
            sum = FIELD[((x + 1) % w) + y * w] + FIELD[((x - 1 + w) % w) + y * w]
                + FIELD[x + ((y + 1) % h) * w] + FIELD[x + ((y - 1 + h) % h) * w];
        } else if (ising_N === 2) {
            sum = FIELD[((x - 1 + w) % w) + ((y - 1 + h) % h) * w] + FIELD[x + ((y - 1 + h) % h) * w]
                + FIELD[((x + 1) % w) + ((y - 1 + h) % h) * w] + FIELD[((x - 1 + w) % w) + y * w]
                + FIELD[((x + 1) % w) + y * w] + FIELD[((x - 1 + w) % w) + ((y + 1) % h) * w]
                + FIELD[x + ((y + 1) % h) * w] + FIELD[((x + 1) % w) + ((y + 1) % h) * w];
        } else {
            for (let dx = -ising_N; dx <= ising_N; dx++) {
                for (let dy = -ising_N; dy <= ising_N; dy++) {
                    if (dx === 0 && dy === 0) continue;
                    const distSq = dx * dx + dy * dy;
                    sum += FIELD[((x + dx + w) % w) + ((y + dy + h) % h) * w] / distSq;
                }
            }
        }

        const deltaE = 2 * spin * (ising_J * sum + ising_B);

        if ((deltaE < 0) || (Math.random() < Math.exp(-deltaE / ising_T)))
        {
            FIELD[i] = -spin;
        }
    }
}

function sim_grayScottStep() {
    const w = GridSizeX;
    const h = GridSizeY;

    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            const i = x + y * w;

            const lapV = lapGrid(V_FIELD, x, y);

            if (Math.abs(lapV) > 1e-6)
            {
                const lapU = lapGrid(U_FIELD, x, y);
                const u = U_FIELD[i];
                const v = V_FIELD[i];
                const uv2 = u * v * v;

                const uNew = u + cahn_dt * (gs_Du * lapU - uv2 + gs_F * (1 - u));
                const vNew = v + cahn_dt * (gs_Dv * lapV + uv2 - (gs_F + gs_k) * v);

                U_NEXT[i] = uNew < 0 ? 0 : (uNew > 1 ? 1 : uNew);
                 V_NEXT[i] = vNew < 0 ? 0 : (vNew > 1 ? 1 : vNew);
             }
         }
     }

     [U_FIELD, U_NEXT] = [U_NEXT, U_FIELD];
     [V_FIELD, V_NEXT] = [V_NEXT, V_FIELD];
 }

 function sim_kobayashiStep() {
     const w = GridSizeX;
     const h = GridSizeY;
     
     // Compute the variational derivative: f'(phi) - epsilon^2 * laplacian(phi) + g(phi)*undercooling
     for (let y = 0; y < h; y++) {
         for (let x = 0; x < w; x++) {
             const i = x + y * w;
             const phi = FIELD[i];
             
             // Double-well potential derivative: f'(phi) = 2 * phi * (1 - phi) * (1 - 2 * phi)
             const dwdphi = 2 * phi * (1 - phi) * (1 - 2 * phi);
             
             // Laplacian term
             const lap = lapGrid(FIELD, x, y);
             
             // Interpolation function g(phi) = phi^2 * (3 - 2*phi)
             // This gives g(0)=0, g(1)=1, smooth transition
             const g_phi = phi * phi * (3 - 2 * phi);
             
             // Variational derivative: f'(phi) - epsilon^2 * laplacian(phi) + g(phi)*undercooling
             const variational = dwdphi - kobayashi_epsilon * kobayashi_epsilon * lap + g_phi * kobayashi_undercooling;
             
             // Time evolution: dphi/dt = -mobility * variational
             const dphi = -kobayashi_mobility * variational;
             
             // Update field
             let newVal = phi + kobayashi_dt * dphi;
             // Keep phi in reasonable bounds [0,1]
             FIELD_NEXT[i] = newVal < 0 ? 0 : (newVal > 1 ? 1 : newVal);
         }
     }
     // Swap buffers
     [FIELD, FIELD_NEXT] = [FIELD_NEXT, FIELD];
}


// DLA helper functions — DO NOT MODIFY (specific to DLA model)
function atGrid_Y(grid, x, y) {
    if (x < 0 || x >= GridSizeX) return 0;
    return grid[x][(y + GridSizeY) % GridSizeY];
}

function isNearCluster(x, y) {
    for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
            if (dx === 0 && dy === 0) continue;
            if (atGrid_Y(FIELD, x + dx, y + dy) === 1) {
                return true;
            }
        }
    }
    return false;
}

function sim_dlaStep() {
    const w = GridSizeX;
    for (let p = 0; p < dla_N; p++) {
        let px = dla_Rightmost + dla_D;
        let py = Math.floor(Math.random() * GridSizeY);

        if (px >= w) continue;

        let steps = 0;
        while (steps < 5000) {
            let near = false;
            for (let dy = -1; dy <= 1 && !near; dy++) {
                for (let dx = -1; dx <= 1 && !near; dx++) {
                    if (dx === 0 && dy === 0) continue;
                    const checkIdx = px + dx + (py + dy) * w;
                    if (checkIdx >= 0 && checkIdx < w * GridSizeY && FIELD[checkIdx] === 1) {
                        near = true;
                        break;
                    }
                }
            }

            if (near && Math.random() < dla_P) {
                if (px >= 0 && px < w) {
                    FIELD[px + py * w] = 1;
                    if (px > dla_Rightmost) dla_Rightmost = px;
                }
                break;
            }

            const dirIdx = Math.floor(Math.random() * 4);
            const newX = px + cartDirections[dirIdx][0];
            const newY = (py + cartDirections[dirIdx][1] + GridSizeY) % GridSizeY;

            if (newX < 0 || newX >= w) break;

            px = newX;
            py = newY;
            steps++;
        }
    }
}


// ============================================================================
// ============================================================================
// ============================================================================
// MODEL INIT FUNCTIONS — EXTENSIBLE
// ============================================================================
// ============================================================================
// ============================================================================
// Each model has an init function that allocates field arrays and sets
// initial conditions. Called when the model is selected or reset.
//
// RULES:
//   - Signature: () => void  (no arguments)
//   - Read GridSizeX, GridSizeY for array dimensions
//   - Use allocIfNeeded(arr, size) instead of new Float32Array(size)
//   - Check `const fresh = !FIELD || FIELD.length !== size` to detect
//     whether a new allocation occurred (grid resize or first call)
//   - Only seed initial values when fresh === true, otherwise reuse
//     existing field data (preserves values across compatible model switches)
//   - Set any model-specific default parameters (e.g., cahn_dt = 1.0)
//   - Do NOT call sim_drawGrid() — the caller handles that
//
// FIELD PRESERVATION on model switch:
//   Cahn ↔ Allen: preserves FIELD (both use φ ∈ [0,1])
//   Cahn/Allen → Ising: converts FIELD [0,1]→[-1,1] before init
//   Ising → Cahn/Allen: converts FIELD [-1,1]→[0,1] before init
//   DLA values (0/1): compatible with both domains, reused as-is
//   All other switches: FIELD is nulled by sim_Reset, re-seeded on init
//
// Available field slots:
//   FIELD, FIELD_NEXT          — single-field models (Cahn, Allen, Ising, DLA)
//   FIELD_MU                   — chemical potential (Cahn-Hilliard)
//   U_FIELD, V_FIELD,          — dual-field models (Gray-Scott)
//   U_NEXT, V_NEXT
//
// If your model needs its own field arrays, declare them in the Field
// Declarations section at the top of this file.
// ============================================================================

function initCahn() {
    const size = GridSizeX * GridSizeY;
    const fresh = !FIELD || FIELD.length !== size;
    FIELD       = allocIfNeeded(FIELD, size);
    FIELD_NEXT  = allocIfNeeded(FIELD_NEXT, size);
    FIELD_MU    = allocIfNeeded(FIELD_MU, size);
    if (fresh) {
        for (let i = 0; i < size; i++) {
            FIELD[i] = 0.5 + (Math.random() - 0.5) * 0.4;
        }
    }
}

function initAllen() {
    const size = GridSizeX * GridSizeY;
    const fresh = !FIELD || FIELD.length !== size;
    FIELD       = allocIfNeeded(FIELD, size);
    FIELD_NEXT  = allocIfNeeded(FIELD_NEXT, size);
    if (fresh) {
        for (let i = 0; i < size; i++) {
            FIELD[i] = 0.5 + (Math.random() - 0.5) * 0.4;
        }
    }
}

function initIsing() {
    const size = GridSizeX * GridSizeY;
    const fresh = !FIELD || FIELD.length !== size;
    FIELD = allocIfNeeded(FIELD, size);
    if (fresh) {
        for (let i = 0; i < size; i++) {
            FIELD[i] = Math.random() < 0.5 ? 1 : -1;
        }
    }
}

function initDLA() {
    const size = GridSizeX * GridSizeY;
    const fresh = !FIELD || FIELD.length !== size;
    FIELD = allocIfNeeded(FIELD, size);
    if (fresh) {
        dla_Rightmost = 0;
        for (let y = 0; y < GridSizeY; y++) {
            FIELD[y * GridSizeX] = 1;
        }
    }
}

function initGray() {
     cahn_dt = 1.0;
     cahn_dx = 1.0;
     const size = GridSizeX * GridSizeY;
     const cx = Math.floor(GridSizeX / 2);
     const cy = Math.floor(GridSizeY / 2);

     U_FIELD = allocIfNeeded(U_FIELD, size);
     V_FIELD = allocIfNeeded(V_FIELD, size);
     U_NEXT  = allocIfNeeded(U_NEXT, size);
     V_NEXT  = allocIfNeeded(V_NEXT, size);

     for (let i = 0; i < size; i++) {
         U_FIELD[i] = 1.0;
         V_FIELD[i] = 0.0;
     }

     const rad = Math.floor(81 / 8);
     for (let x = 0; x < GridSizeX; x++) {
         for (let y = 0; y < GridSizeY; y++) {
             if ((x-cx)*(x-cx) + (y-cy)*(y-cy) <= rad*rad) {
                 const i = x + y * GridSizeX;
                 U_FIELD[i] = 0.5 + 0.01 * (Math.random() - 0.5);
                 V_FIELD[i] = 0.25 + 0.01 * (Math.random() - 0.5);
             }
         }
     }
 }

 function initKobayashi() {
     const size = GridSizeX * GridSizeY;
     const fresh = !FIELD || FIELD.length !== size;
     FIELD       = allocIfNeeded(FIELD, size);
     FIELD_NEXT  = allocIfNeeded(FIELD_NEXT, size);
     if (fresh) {
         // Initialize with small solid nucleus in liquid matrix
         for (let i = 0; i < size; i++) {
             FIELD[i] = 0.0;  // Start with all liquid
         }
         // Add a small solid seed at center
         const cx = Math.floor(GridSizeX / 2);
         const cy = Math.floor(GridSizeY / 2);
         const radius = 3;
         for (let y = Math.max(0, cy - radius); y < Math.min(GridSizeY, cy + radius); y++) {
             for (let x = Math.max(0, cx - radius); x < Math.min(GridSizeX, cx + radius); x++) {
                 if ((x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius*radius) {
                     FIELD[x + y * GridSizeX] = 1.0;  // Solid nucleus
                 }
             }
         }
     }
 }


// ============================================================================
// ============================================================================
// ============================================================================
// MODEL DRAW FUNCTIONS — EXTENSIBLE
// ============================================================================
// ============================================================================
// ============================================================================
// Each model has a draw function that maps field values to pixel colors.
//
// RULES:
//   - Signature: () => void  (no arguments)
//   - Must create/reuse simImageData if needed (check dimensions)
//   - Map field values to [0, 255] palette index
//   - Call drawScaleUp() at the end to render to canvas
//   - Do NOT use if/else for model selection — each model gets its own
//     draw function registered in modelConfigs
//
// REUSABLE DRAW FUNCTIONS:
//   drawField01()     — maps FIELD[i] ∈ [0,1] to palette. Reuse for any
//                       model whose field is normalized to [0,1].
//                       Used by: Cahn-Hilliard, Allen-Cahn, DLA.
//   drawFieldIsing()  — maps FIELD[i] ∈ {-1,+1} to palette via (x+1)/2.
//                       Reuse for any model with binary ±1 states.
//   drawFieldGray()   — maps V_FIELD[i] / 0.25 to palette.
//                       Gray-Scott specific.
//
// To add a new draw function, write it here and reference it in modelConfigs.
// ============================================================================

function drawField01() {
    if (!simImageData || simImageData.width !== GridSizeX || simImageData.height !== GridSizeY) {
        simImageData = simCTX.createImageData(GridSizeX, GridSizeY);
    }

    const data = simImageData.data;
    const palette = colorPalette;
    const size = GridSizeX * GridSizeY;

    for (let i = 0; i < size; i++) {
        const val = Math.floor(Math.max(0, Math.min(1, FIELD[i])) * 255);
        const c = val * 3;
        data[i * 4 + 0] = palette[c + 0];
        data[i * 4 + 1] = palette[c + 1];
        data[i * 4 + 2] = palette[c + 2];
        data[i * 4 + 3] = 255;
    }

    drawScaleUp();
}

function drawFieldIsing() {
    if (!simImageData || simImageData.width !== GridSizeX || simImageData.height !== GridSizeY) {
        simImageData = simCTX.createImageData(GridSizeX, GridSizeY);
    }

    const data = simImageData.data;
    const palette = colorPalette;
    const size = GridSizeX * GridSizeY;

    for (let i = 0; i < size; i++) {
        let val = (FIELD[i] + 1) / 2;
        val = Math.max(0, Math.min(1, val));
        const c = Math.floor(val * 255) * 3;
        data[i * 4 + 0] = palette[c + 0];
        data[i * 4 + 1] = palette[c + 1];
        data[i * 4 + 2] = palette[c + 2];
        data[i * 4 + 3] = 255;
    }

    drawScaleUp();
}

function drawFieldGray() {
    if (!simImageData || simImageData.width !== GridSizeX || simImageData.height !== GridSizeY) {
        simImageData = simCTX.createImageData(GridSizeX, GridSizeY);
    }

    const data = simImageData.data;
    const palette = colorPalette;
    const size = GridSizeX * GridSizeY;

    for (let i = 0; i < size; i++) {
        let val = V_FIELD[i] / 0.25;
        val = Math.max(0, Math.min(1, val));
        const c = Math.floor(val * 255) * 3;
        data[i * 4 + 0] = palette[c + 0];
        data[i * 4 + 1] = palette[c + 1];
        data[i * 4 + 2] = palette[c + 2];
        data[i * 4 + 3] = 255;
    }

    drawScaleUp();
}

// CRITICAL: shared by all draw functions. DO NOT MODIFY.
function drawScaleUp() {
    const offCanvas = document.createElement('canvas');
    offCanvas.width = GridSizeX;
    offCanvas.height = GridSizeY;
    const offCtx = offCanvas.getContext('2d');
    offCtx.putImageData(simImageData, 0, 0);

    simCTX.imageSmoothingEnabled = false;
    simCTX.drawImage(offCanvas, 0, 0, canvas.width, canvas.height);
}


// ============================================================================
// FIELD CONVERSION — DO NOT MODIFY
// ============================================================================
// Used when switching between phase-field models (Cahn/Allen) and Ising.
// DLA values (0/1) are already compatible with both domains, so no
// conversion needed — just reuse via allocIfNeeded.
// ============================================================================

function isPhaseModel(m) { return m === 'cahnhilliard' || m === 'allencahn' || m === 'kobayashi'; }
function isIsingModel(m) { return m === 'ising'; }

function convertFieldToIsing() {
    const size = GridSizeX * GridSizeY;
    for (let i = 0; i < size; i++) {
        FIELD[i] = FIELD[i] * 2 - 1;
    }
}

function convertFieldFromIsing() {
    const size = GridSizeX * GridSizeY;
    for (let i = 0; i < size; i++) {
        FIELD[i] = (FIELD[i] + 1) / 2;
    }
}


// ============================================================================
// ============================================================================
// ============================================================================
// MODEL CONFIG — DO NOT MODIFY STRUCTURE — EXTENSIBLE ENTRIES
// ============================================================================
// ============================================================================
// ============================================================================
// CRITICAL: This is the central dispatch table. The simulation engine uses
// this to find the right init/step/draw functions for each model.
//
// STRUCTURE — DO NOT MODIFY:
//   modelConfigs: { [modelKey]: { init: Function, step: Function, draw: Function } }
//   sim_initFn, sim_stepFn, sim_drawFn: set from modelConfigs on model select
//
// TO ADD A NEW MODEL — ADD ONE ENTRY:
//   yourModelKey: { init: initYourModel, step: sim_yourModelStep, draw: drawYourModel }
//
// The modelKey MUST match:
//   1. The data-model attribute in index.html dropdown item
//   2. The key used in paramConfigs (see Section: Param Configs)
//
// draw reuse guide:
//   - Field ∈ [0,1] single array?  use drawField01
//   - Field ∈ {-1,+1}?            use drawFieldIsing
//   - Dual-field (U,V)?           write a new draw function
// ============================================================================

const modelConfigs = {
     cahnhilliard: { init: initCahn,  step: sim_cahnHilliardStep, draw: drawField01 },
     allencahn:    { init: initAllen, step: sim_allenCahnStep,    draw: drawField01 },
     ising:        { init: initIsing, step: sim_isingStep,        draw: drawFieldIsing },
     dla:          { init: initDLA,   step: sim_dlaStep,          draw: drawField01 },
     gray:         { init: initGray,  step: sim_grayScottStep,    draw: drawFieldGray },
     kobayashi:    { init: initKobayashi, step: sim_kobayashiStep, draw: drawField01 },
 };

// CRITICAL: global function references. Set by model selection. DO NOT MODIFY.
let sim_initFn = modelConfigs[simModel].init;
let sim_stepFn = modelConfigs[simModel].step;
let sim_drawFn = modelConfigs[simModel].draw;


// ============================================================================
// ============================================================================
// ============================================================================
// SIMULATION RUN LOOP — DO NOT MODIFY
// ============================================================================
// ============================================================================
// ============================================================================
// CRITICAL: sim_Run is the animation loop. It calls sim_stepFn() then
// sim_drawFn() via requestAnimationFrame. Do not add if/else or model
// checks here — all model dispatch goes through the function references.
// ============================================================================

function sim_Run(timestamp) {
    sim_stepFn();
    sim_drawFn();

    if (simRunning) {
        simRequest = requestAnimationFrame(sim_Run);
    }
}


// ============================================================================
// SIMULATION CONTROLS — DO NOT MODIFY
// ============================================================================

function sim_Start() {
    if (!simRunning) {
        simRunning = true;
        simRequest = requestAnimationFrame(sim_Run);
        document.getElementById('startBtn').disabled = true;
        document.getElementById('stopBtn').disabled  = false;
        document.getElementById('resetBtn').disabled = true;
    }
}

function sim_Stop() {
    if (simRunning) {
        cancelAnimationFrame(simRequest);
        simRunning = false;
        document.getElementById('startBtn').disabled = false;
        document.getElementById('stopBtn').disabled  = true;
        document.getElementById('resetBtn').disabled = false;
    }
}

function sim_Reset() {
    sim_Stop();
    FIELD = null; FIELD_NEXT = null; FIELD_MU = null;
    U_FIELD = null; V_FIELD = null; U_NEXT = null; V_NEXT = null;
    sim_initModelGrid();
    sim_drawGrid();
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled  = true;
    document.getElementById('resetBtn').disabled = false;
}

function sim_UpdateGridSize(event) {
    sim_Stop();

    GridSizeX = parseInt(event.target.value);
    GridSizeY = sim_GetSizeY(GridSizeX);

    document.getElementById('GridSizeLabel').innerText = `Grid Size: ${GridSizeX} x ${GridSizeY}`;

    sim_initModelGrid();
    sim_drawGrid();

    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled = true;
    document.getElementById('resetBtn').disabled = false;
}


// ============================================================================
// SLIDER UI SYSTEM — DO NOT MODIFY
// ============================================================================
// sim_createSliderControl builds a slider from a config object.
// sim_SetInputParmeters reads paramConfigs[simModel] to populate sliders.
// ============================================================================

function sim_createSliderControl(paramDiv, cfg)
{
    const id = cfg.id;
    const min = cfg.min ?? 0;
    const max = cfg.max ?? 1;
    const step = cfg.step ?? 0.01;
    const toFixed = cfg.toFixed ?? 2;
    const getValue = cfg.getValue;
    const onChange = cfg.onChange ?? (() => {});

    const row         = document.createElement('div');
    row.className     = 'parameter-row';

    const labelDiv    = document.createElement('div');
    labelDiv.className = 'parameter-label';

    const label       = document.createElement('span');
    label.textContent = cfg.label;

    const valueSpan   = document.createElement('span');
    valueSpan.id      = `${id}Value`;
    valueSpan.className = 'parameter-value';
    const value = getValue();
    valueSpan.textContent = (toFixed !== undefined) ? value.toFixed(toFixed) : value.toExponential(6);

    labelDiv.appendChild(label);
    labelDiv.appendChild(valueSpan);

    const slider      = document.createElement('input');
    slider.id         = id;
    slider.type       = 'range';
    slider.min        = min;
    slider.max        = max;
    slider.step       = step;
    slider.value      = value;
    slider.className = 'sidebar-slider';

    slider.addEventListener('input', (e) => {
        const newValue = parseFloat(e.target.value);
        onChange(newValue);
        valueSpan.textContent = (toFixed !== undefined) ? newValue.toFixed(toFixed) : newValue.toExponential(6);
    });

    row.appendChild(labelDiv);
    row.appendChild(slider);
    paramDiv.appendChild(row);
}


// ============================================================================
// ============================================================================
// ============================================================================
// PARAM CONFIGS — DO NOT MODIFY STRUCTURE — EXTENSIBLE ENTRIES
// ============================================================================
// ============================================================================
// ============================================================================
// Each model has an array of slider config objects. These are rendered by
// sim_SetInputParmeters when a model is selected.
//
// The key for each entry MUST match the corresponding key in modelConfigs.
// The cellSize entry is shared by all models — do not remove it.
//
// TO ADD A NEW MODEL — ADD ONE ENTRY:
//   yourModelKey: [
//       { id, label, min, max, step, toFixed, getValue, onChange },
//       ...
//   ]
//
// Config object shape:
//   id        — unique DOM id for the slider
//   label     — display label shown above the slider
//   min/max   — slider range
//   step      — slider step size
//   toFixed   — decimal places for value display (use 0 for integers)
//   getValue  — () => number, reads the current parameter value
//   onChange  — (newValue) => void, sets the parameter value
// ============================================================================

const paramConfigs = {
    cellSize: {
        id: 'cellSize',
        label: 'Cell Size',
        min: 1,
        max: 20,
        step: 1,
        toFixed: 0,
        getValue: () => cellSize,
        onChange: (val) => { cellSize = val; sim_resizeCanvas(); }
    },
    cahnhilliard: [
        { id: 'dxLabel', label: 'Grid Spacing (Δx)', min: 0.0, max: 5, step: 0.1, toFixed: 2, getValue: () => cahn_dx, onChange: (val) => cahn_dx = val },
        { id: 'dt',      label: 'Time Step (Δt)',    min: 0.0, max: 1, step: 0.1, toFixed: 2, getValue: () => cahn_dt, onChange: (val) => cahn_dt = val },

        { id: 'M',       label: 'Mobility (M)',             min: 0.0, max: 1.0, step: 0.1, toFixed: 2, getValue: () => cahn_M, onChange: (val) => cahn_M = val },
        { id: 'k',       label: 'Interfacial Energy (κ)',   min: 0.0, max: 1.0, step: 0.1, toFixed: 2, getValue: () => cahn_K, onChange: (val) => cahn_K = val },
        { id: 'W',       label: 'Free Energy Coeff (A)',    min: 0.0, max: 1.0, step: 0.1, toFixed: 2, getValue: () => cahn_W, onChange: (val) => cahn_W = val }
    ],
    ising: [
        { id: 'T', label: 'Temperature (T)',        min: 0, max: 5, step: 0.05, toFixed: 2, getValue: () => ising_T, onChange: (val) => ising_T = val },
        { id: 'H', label: 'Magnetic Field (H)',     min:-2, max: 2, step: 0.05, toFixed: 2, getValue: () => ising_B, onChange: (val) => ising_B = val },
        { id: 'J', label: 'Coupling constant (J)',  min: 0, max: 5, step: 0.05, toFixed: 2, getValue: () => ising_J, onChange: (val) => ising_J = val },
        { id: 'N', label: 'Neighborhood (N)',       min: 1, max: 5, step: 1,    toFixed: 2, getValue: () => ising_N, onChange: (val) => ising_N = val }
    ],
    dla: [
        { id: 'P',  label: 'Particles per Step',    min: 0, max: 1000, step: 50,   toFixed: 0, getValue: () => dla_N, onChange: (val) => dla_N = val },
        { id: 'P2', label: 'Sticking Probability',  min: 0, max: 1.0,  step: 0.01, toFixed: 2, getValue: () => dla_P, onChange: (val) => dla_P = val },
        { id: 'D',  label: 'Launch Distance',       min: 0, max: 20,   step: 1,    toFixed: 0, getValue: () => dla_D, onChange: (val) => dla_D = val }
    ],
    gray: [
        { id: 'dxLabel', label: 'Grid Spacing (Δx)', min: 0.0, max: 5, step: 0.01, toFixed: 2, getValue: () => cahn_dx, onChange: (val) => cahn_dx = val },
        { id: 'dt',      label: 'Time Step (Δt)',    min: 0.0, max: 1, step: 0.01, toFixed: 2, getValue: () => cahn_dt, onChange: (val) => cahn_dt = val },

        { id: 'Du', label: 'Diffusion U (D_u)', min: 0.0, max: 1.0, step: 0.01, toFixed: 2, getValue: () => gs_Du, onChange: (val) => gs_Du = val },
        { id: 'Dv', label: 'Diffusion V (D_v)', min: 0.0, max: 1.0, step: 0.01, toFixed: 2, getValue: () => gs_Dv, onChange: (val) => gs_Dv = val },
        { id: 'F',  label: 'Feed Rate (F)',     min: 0.0, max: 1.0, step: 0.01, toFixed: 2, getValue: () => gs_F, onChange: (val) => gs_F = val },
        { id: 'k',  label: 'Kill Rate (k)',     min: 0.0, max: 1.0, step: 0.01, toFixed: 2, getValue: () => gs_k, onChange: (val) => gs_k = val }
    ],
     allencahn: [
         { id: 'ac_dx', label: 'Grid Spacing (Δx)', min: 0.0, max: 5, step: 0.1, toFixed: 2, getValue: () => cahn_dx, onChange: (val) => cahn_dx = val },
         { id: 'ac_dt', label: 'Time Step (Δt)',    min: 0.0, max: 1, step: 0.1, toFixed: 2, getValue: () => cahn_dt, onChange: (val) => cahn_dt = val },

         { id: 'ac_M', label: 'Mobility (M)',              min: 0.0, max: 2.0, step: 0.1, toFixed: 2, getValue: () => ac_M, onChange: (val) => ac_M = val },
         { id: 'ac_K', label: 'Interfacial Energy (κ)',    min: 0.0, max: 2.0, step: 0.1, toFixed: 2, getValue: () => ac_K, onChange: (val) => ac_K = val },
         { id: 'ac_W', label: 'Free Energy Coeff (W)',     min: 0.0, max: 2.0, step: 0.1, toFixed: 2, getValue: () => ac_W, onChange: (val) => ac_W = val },
         { id: 'ac_F', label: 'Driving Force (F)',         min:-0.1, max: 0.1, step: 0.01, toFixed: 2, getValue: () => ac_F, onChange: (val) => ac_F = val }
     ],
     kobayashi: [
         { id: 'kob_dt', label: 'Time Step (Δt)', min: 0.0, max: 0.5, step: 0.01, toFixed: 3, getValue: () => kobayashi_dt, onChange: (val) => kobayashi_dt = val },
         { id: 'kob_epsilon', label: 'Interface Thickness (ε)', min: 0.1, max: 2.0, step: 0.1, toFixed: 2, getValue: () => kobayashi_epsilon, onChange: (val) => kobayashi_epsilon = val },
         { id: 'kob_M', label: 'Mobility (M)', min: 0.0, max: 2.0, step: 0.1, toFixed: 2, getValue: () => kobayashi_mobility, onChange: (val) => kobayashi_mobility = val },
         { id: 'kob_dT', label: 'Undercooling (ΔT)', min: -1.0, max: 1.0, step: 0.05, toFixed: 2, getValue: () => kobayashi_undercooling, onChange: (val) => kobayashi_undercooling = val }

     ]
 };


// ============================================================================
// PARAMETER SETUP — DO NOT MODIFY
// ============================================================================
// CRITICAL: sim_SetInputParmeters is the slider dispatcher. It reads
// paramConfigs[simModel] and populates the sidebar. Do not add if/else here.
// ============================================================================

function sim_SetInputParmeters() {
    const paramDiv = document.getElementById('parameters');
    paramDiv.innerHTML = '';

    const config = paramConfigs[simModel] || [];

    sim_createSliderControl(paramDiv, paramConfigs.cellSize);

    config.forEach(cfg => sim_createSliderControl(paramDiv, cfg));
}


// ============================================================================
// ============================================================================
// ============================================================================
// EVENT WIRING — DO NOT MODIFY
// ============================================================================
// ============================================================================
// ============================================================================
// CRITICAL: window.onload sets up all event listeners and initializes the
// simulation. The dropdown handler sets simModel and refreshes the global
// function references (sim_initFn, sim_stepFn, sim_drawFn).
// Do not modify this section.
// ============================================================================

window.onload = () => {
    const menuToggle      = document.getElementById('menuToggle');
    const dropdownContent = document.getElementById('dropdownContent');

    menuToggle.addEventListener('click', function(event) {
        dropdownContent.classList.toggle('show');
        menuToggle.classList.toggle('active');
        event.stopPropagation();
    });

    dropdownContent.querySelectorAll('.dropdown-item').forEach(item => {
        item.addEventListener('click', function(event) {
            event.preventDefault();
            stopWelcomeAnimation();
            const selectedModel = this.getAttribute('data-model');
            const prevModel = simModel;

            simModel = selectedModel;

            sim_initFn = modelConfigs[simModel].init;
            sim_stepFn = modelConfigs[simModel].step;
            sim_drawFn = modelConfigs[simModel].draw;

            menuToggle.childNodes[0].textContent = this.textContent + ' ';

            sim_Stop();

            // Convert FIELD before init when switching between compatible domains
            if (simModel === 'dla') {
                FIELD = null;
            } else if (isPhaseModel(prevModel) && isIsingModel(simModel)) {
                const size = GridSizeX * GridSizeY;
                FIELD = allocIfNeeded(FIELD, size);
                convertFieldToIsing();
            } else if (isIsingModel(prevModel) && isPhaseModel(simModel)) {
                const size = GridSizeX * GridSizeY;
                FIELD = allocIfNeeded(FIELD, size);
                convertFieldFromIsing();
            }

            sim_initFn();
            sim_SetInputParmeters();
            sim_drawGrid();

            dropdownContent.classList.remove('show');
            menuToggle.classList.remove('active');
        });
    });


    window.addEventListener('click', function(event) {
        if (!event.target.closest('.dropdown-menu-container')) {
            if (dropdownContent.classList.contains('show')) {
                dropdownContent.classList.remove('show');
                menuToggle.classList.remove('active');
            }
        }
    });

    document.getElementById('startBtn').addEventListener('click', sim_Start);
    document.getElementById('stopBtn').addEventListener('click', sim_Stop);
    document.getElementById('resetBtn').addEventListener('click', sim_Reset);
    const gridSizeInput = document.getElementById('GridSizeInput');
    if (gridSizeInput) {
        gridSizeInput.addEventListener('input', sim_UpdateGridSize);
    }

    // Cell size slider
    const cellSizeSlider = document.getElementById('cellSizeSlider');
    const cellSizeValue = document.getElementById('cellSizeValue');
    if (cellSizeSlider) {
        cellSize = parseInt(cellSizeSlider.value);
        cellSizeSlider.addEventListener('input', (e) => {
            cellSize = parseInt(e.target.value);
            if (cellSizeValue) {
                cellSizeValue.textContent = cellSize;
            }
            sim_resizeCanvas();
        });
    }

    window.addEventListener('resize', sim_resizeCanvas);

    initColorPalette();

    // Defer initialization to ensure DOM is ready
    setTimeout(() => {
        sim_resizeCanvas();
        if (USE_MAXIMALIST_WELCOME) {
            initWelcomeParticles();
        }
        displayWelcomeCanvas();
        sim_Stop();
    }, 100);
};

function stopWelcomeAnimation() {
    welcomeRunning = false;
    if (welcomeAnimId) {
        cancelAnimationFrame(welcomeAnimId);
        welcomeAnimId = null;
    }
}

function restartWelcomeAnimation() {
    if (welcomeAnimId) {
        cancelAnimationFrame(welcomeAnimId);
    }
    welcomeAnimFrame = 0;
    welcomeRunning = true;
    animateWelcome();
}
