let currentModel = 'cahn-hilliard';
let ctx = document.getElementById('canvas').getContext('2d');

let phi;
let isingGrid;
let simulationInterval;

let deltaX   = 1.0;
let gridSize = 128;
let running  = false;

let ch_timeStep = 0.02;  
let ch_mobility = 0.3;
let ch_kappa    = 0.3;
let ch_Wcoeff   = 1.0;

let ising_T = 2.269; 
let ising_B    = 0.0; 
let ising_J    = 1.0; 
const is_Kb = 1.0;

function sim_resizeCanvas()
{
    const canvas  = document.getElementById('canvas');
    canvas.width  = window.innerWidth - 320;
    canvas.height = window.innerHeight;

    const minCanvasDim = Math.min(canvas.width, canvas.height);
    const maxGridSize  = Math.floor(minCanvasDim / 3); 
    
    const gridSizeInput = document.getElementById('gridSizeInput');
    gridSizeInput.setAttribute('max', maxGridSize);

    if (gridSize > maxGridSize) {
        gridSize = maxGridSize;
        gridSizeInput.value = gridSize;
    }
    
    document.getElementById('gridSizeLabel').innerText = `Size: ${gridSize} x ${gridSize}`;
    sim_initModelGrid();
    sim_drawGrid();
}

function sim_initModelGrid()
{
    sim_Stop();

    if (currentModel === 'cahn-hilliard')
    {
        phi = Array.from({ length: gridSize }, () => 
            Array.from({ length: gridSize }, () => 0.5 + (Math.random() - 0.5) * 0.4)
        );

    } else if (currentModel === 'ising')
    {
        isingGrid = Array.from({ length: gridSize }, () =>
            Array.from({ length: gridSize }, () => (Math.random() < 0.5 ? 1 : -1))
        );
    }
}

function hexToRgb(hex) {
    const bigint = parseInt(hex.slice(1), 16);
    const r = (bigint >> 16) & 255;
    const g = (bigint >> 8) & 255;
    const b = bigint & 255;
    return [r, g, b];
}

// Generated
function sim_drawGrid()
{
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    const renderCellSize = ctx.canvas.width / gridSize;

    // Define colors for both models
    const colorHexA = "#0D5EA6";
    const colorHexB = "#EAA64D"; 

    const rgbA = hexToRgb(colorHexA);
    const rgbB = hexToRgb(colorHexB);

    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {
            let normalizedValue;

            if (currentModel === 'cahn-hilliard') {
                normalizedValue = phi[x][y];
            } else if (currentModel === 'ising') {
                normalizedValue = (isingGrid[x][y] + 1) / 2; 
            }

            const red   = Math.floor(rgbA[0] + (rgbB[0] - rgbA[0]) * normalizedValue);
            const green = Math.floor(rgbA[1] + (rgbB[1] - rgbA[1]) * normalizedValue);
            const blue  = Math.floor(rgbA[2] + (rgbB[2] - rgbA[2]) * normalizedValue);

            ctx.fillStyle = `rgb(${red}, ${green}, ${blue})`;
            ctx.fillRect(x * renderCellSize, y * renderCellSize, renderCellSize, renderCellSize);
        }
    }
}

function atGrid(grid, x, y)
{
    return grid[(x + gridSize) % gridSize][(y + gridSize) % gridSize];
}

function lapGrid(field, x, y)
{
    return (atGrid(field, x + 1, y) + atGrid(field, x - 1, y) +
            atGrid(field, x, y + 1) + atGrid(field, x, y - 1) - 4 * atGrid(field, x, y)) / (deltaX * deltaX);
}

function sim_cahnHilliardStep()
{
    const mu     = Array.from({ length: gridSize }, () => Array(gridSize).fill(0));
    const newPhi = Array.from({ length: gridSize }, () => Array(gridSize).fill(0));

    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {
            const phiVal = atGrid(phi, x, y);
            const f_prime_phi = 2 * ch_Wcoeff * phiVal * (1 - phiVal) * (1 - 2 * phiVal);
            mu[x][y] = f_prime_phi - ch_kappa * lapGrid(phi, x, y);
        }
    }

    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {
            const lapMu = lapGrid(mu, x, y);
            newPhi[x][y] = atGrid(phi, x, y) + ch_timeStep * ch_mobility * lapMu;
            newPhi[x][y] = Math.max(0, Math.min(1, newPhi[x][y]));
        }
    }

    phi = newPhi;
}

function sim_isingStep()
{
    for (let i = 0; i < gridSize * gridSize; i++) {

        const x = Math.floor(Math.random() * gridSize);
        const y = Math.floor(Math.random() * gridSize);

        const currentSpin = atGrid(isingGrid, x, y);

        // Calculate sum of neighbor spins
        const sumNeighbors = atGrid(isingGrid, x + 1, y) +
                             atGrid(isingGrid, x - 1, y) +
                             atGrid(isingGrid, x, y + 1) +
                             atGrid(isingGrid, x, y - 1);

        const deltaE = 2 * ising_J * currentSpin * sumNeighbors + 2 * ising_B * currentSpin;

        if (deltaE < 0) {
            isingGrid[x][y] *= -1;
        } else {
            const p = Math.exp(-deltaE / (is_Kb * ising_T));
            if (Math.random() < p) {
                isingGrid[x][y] *= -1;
            }
        }
    }
}

function sim_Run(){
    if (currentModel === 'cahn-hilliard') {
        sim_cahnHilliardStep();
    } else if (currentModel === 'ising') {
        sim_isingStep();
    }
    sim_drawGrid();
}

function sim_Start() {
    if (!running) {
        simulationInterval = setInterval(sim_Run, 10);
        running = true;
        document.getElementById('startBtn').disabled = true;
        document.getElementById('stopBtn').disabled = false;
        document.getElementById('resetBtn').disabled = true;
    }
}

function sim_Stop() {
    clearInterval(simulationInterval);
    running = false;
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled = true;
    document.getElementById('resetBtn').disabled = false;
}

function sim_Reset() {
    sim_Stop();
    sim_initModelGrid();
    sim_drawGrid();
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled = true;
    document.getElementById('resetBtn').disabled = false;
}

function sim_UpdateGridSize(event) {
    sim_Stop();
    gridSize = parseInt(event.target.value);
    document.getElementById('gridSizeLabel').innerText = `Size: ${gridSize} x ${gridSize}`;
    sim_initModelGrid();
    sim_drawGrid();
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled = true;
    document.getElementById('resetBtn').disabled = false;
}

function sim_createSliderControl(paramDiv, { id, labelPrefix, min, max, step, value, toFixed, onChange }) 
{
    const label         = document.createElement('label');
    label.id            = `${id}Label`;
    label.className     = 'sidebar-label';
    label.textContent   = `${labelPrefix} = ${toFixed ? value.toFixed(toFixed) : value.toExponential(6)}`;

    const slider        = document.createElement('input');
    slider.id           = id;
    slider.type         = 'range';
    slider.min          = min;
    slider.max          = max;
    slider.step         = step;
    slider.value        = value;
    slider.className    = 'sidebar-slider';

    slider.addEventListener('input', (e) => {
        const newValue = parseFloat(e.target.value);
        onChange(newValue);
        label.textContent = `${labelPrefix} = ${toFixed ? newValue.toFixed(toFixed) : newValue.toExponential(6)}`;
        // if (running) sim_Stop();
    });

    paramDiv.appendChild(label);
    paramDiv.appendChild(slider);
}

function sim_SetInputParmeters() {
    const paramDiv = document.getElementById('parameters');
    paramDiv.innerHTML = '';

    if (currentModel === 'cahn-hilliard') {
        
        sim_createSliderControl(paramDiv, {
            id: 'dxLabel',
            labelPrefix: 'Grid Spacing (Δx)',
            min: 0.1,
            max: 5,
            step: 0.01,
            value: deltaX,
            toFixed: 2,
            onChange: (val) => { deltaX = val; }
        });

        sim_createSliderControl(paramDiv, {
            id: 'dt',
            labelPrefix: 'Time sim_Run (Δt)',
            min: 0.001,
            max: 1,
            step: 0.001,
            value: ch_timeStep,
            toFixed: 4, // Exponential display
            onChange: (val) => { ch_timeStep = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'M',
            labelPrefix: 'Mobility (M)',
            min: 0.01,
            max: 1.0,
            step: 0.005,
            value: ch_mobility,
            toFixed: 2,
            onChange: (val) => { ch_mobility = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'k',
            labelPrefix: 'Interfacial Energy (κ)',
            min: 0.01,
            max: 1.0,
            step: 0.005,
            value: ch_kappa,
            toFixed: 2,
            onChange: (val) => { ch_kappa = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'W',
            labelPrefix: 'Free Energy Coeff (A)',
            min: 0.1,
            max: 2.0,
            step: 0.1,
            value: ch_Wcoeff,
            toFixed: 2,
            onChange: (val) => { ch_Wcoeff = val; }
        });

    } else if (currentModel === 'ising') {

        sim_createSliderControl(paramDiv, {
            id: 'T',
            labelPrefix: 'Temperature (T)',
            min: 0,
            max: 5,
            step: 0.01,
            value: ising_T,
            toFixed: 2,
            onChange: (val) => { ising_T = val; }
        });

        sim_createSliderControl(paramDiv, {
            id: 'H',
            labelPrefix: 'Magnetic Field (H)',
            min: -2,
            max: 2,
            step: 0.01,
            value: ising_B,
            toFixed: 2,
            onChange: (val) => { ising_B = val; }
        });

        sim_createSliderControl(paramDiv, {
            id: 'J',
            labelPrefix: 'Coupling constant (J)',
            min: 0,
            max: 5,
            step: 0.01,
            value: ising_J,
            toFixed: 2,
            onChange: (val) => { ising_J = val; }
        });
    }
}


window.onload = () => {

    const menuToggle      = document.getElementById('menuToggle');
    const dropdownContent = document.getElementById('dropdownContent');

    menuToggle.addEventListener('click', function(event) {
        dropdownContent.classList.toggle('show');
        menuToggle.classList.toggle('active');
        event.stopPropagation();
    });

    // Handle dropdown item clicks
    dropdownContent.querySelectorAll('.dropdown-item').forEach(item => {
        item.addEventListener('click', function(event) {
            event.preventDefault();
            const selectedModel = this.getAttribute('data-model');
            // if (selectedModel !== currentModel) {
                currentModel = selectedModel;
                // Update header text to reflect chosen model
                menuToggle.firstChild.textContent = this.textContent + ' '; // Add space back for icon
                
                sim_Stop();
                sim_initModelGrid();
                sim_SetInputParmeters();
                sim_drawGrid();
            // }
            dropdownContent.classList.remove('show'); // Close dropdown
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
    document.getElementById('gridSizeInput').addEventListener('input', sim_UpdateGridSize);
    window.addEventListener('resize', sim_resizeCanvas);

    sim_resizeCanvas();
    // sim_SetInputParmeters();
    sim_Stop();
};