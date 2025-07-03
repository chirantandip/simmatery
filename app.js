let simCTX  = document.getElementById('canvas').getContext('2d');
let simModel;
let simRequest;

const simCahn  = 'cahn-hilliard';
const simIsing = 'ising';

let FIELD;
let FIELD_NEXT;
let FIELD_MU;

let gridSize = 128;
let running  = false;

let cahn_dx = 1.0;
let cahn_dt = 0.25;  
let cahn_M  = 0.25; /// Mobility
let cahn_K  = 0.25; /// Kappa
let cahn_W  = 0.99; /// Double well energy

let ising_T = 2.0;  /// Temperature
let ising_B = 0.0;  /// Magnetif field
let ising_J = 1.2;  /// Coupling constant
let ising_N = 1;    /// neighborhood

function sim_resizeCanvas() {
    const canvas  = document.getElementById('canvas');
    canvas.height = window.innerHeight;
    canvas.width  = window.innerWidth - 400;

    const maxGridSize   = Math.floor(Math.min(canvas.width, canvas.height) / 3); 
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

function sim_initModelGrid() {
    sim_Stop();
    if (simModel === simCahn) {
        FIELD       = Array.from({ length: gridSize }, () => Array.from({ length: gridSize }, () => 0.5 + (Math.random() - 0.5) * 0.4));
        FIELD_NEXT  = Array.from({ length: gridSize }, () => Array(gridSize).fill(0));
        FIELD_MU    = Array.from({ length: gridSize }, () => Array(gridSize).fill(0));

    } else if (simModel === simIsing) {
        FIELD       = Array.from({ length: gridSize }, () => Array.from({ length: gridSize }, () => (Math.random() < 0.5 ? 1 : -1)));
    }
}

function sim_drawGrid() {
    simCTX.clearRect(0, 0, simCTX.canvas.width, simCTX.canvas.height);

    const cellSize = simCTX.canvas.width / gridSize;

    const rgbA = [13, 94, 166];
    const rgbB = [234,166, 77];
    let val;

    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {

            if (simModel === simCahn) {
                val = FIELD[x][y];
            } else if (simModel === simIsing) {
                val = (FIELD[x][y] + 1) / 2; 
            }

            const red   = Math.floor(rgbA[0] + (rgbB[0] - rgbA[0]) * val);
            const green = Math.floor(rgbA[1] + (rgbB[1] - rgbA[1]) * val);
            const blue  = Math.floor(rgbA[2] + (rgbB[2] - rgbA[2]) * val);

            simCTX.fillStyle = `rgb(${red}, ${green}, ${blue})`;
            simCTX.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
        }
    }
}

function atGrid(grid, x, y) {
    return grid[(x + gridSize) % gridSize][(y + gridSize) % gridSize];
}

function lapGrid(field, x, y) {
    return (atGrid(field, x + 1, y) + atGrid(field, x - 1, y) +
            atGrid(field, x, y + 1) + atGrid(field, x, y - 1) - 4 * atGrid(field, x, y)) / (cahn_dx * cahn_dx);
}

function sim_cahnHilliardStep() {
    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {
            const phi       = atGrid(FIELD, x, y);
            const dwdphi    = 2 * cahn_W * phi * (1 - phi) * (1 - 2 * phi);
            FIELD_MU[x][y]  = dwdphi - cahn_K * lapGrid(FIELD, x, y);
        }
    }
    for (let x = 0; x < gridSize; x++) {
        for (let y = 0; y < gridSize; y++) {
            FIELD_NEXT[x][y] = atGrid(FIELD, x, y) + cahn_dt * cahn_M * lapGrid(FIELD_MU, x, y);
            FIELD_NEXT[x][y] = Math.max(0, Math.min(1, FIELD_NEXT[x][y]));
        }
    }

    const tmp = FIELD;
    FIELD = FIELD_NEXT;
    FIELD_NEXT = tmp;
}

function sim_isingStep() {
    for (let i = 0; i < gridSize * gridSize; i++) {
        const x = Math.floor(Math.random() * gridSize);
        const y = Math.floor(Math.random() * gridSize);

        const spin = atGrid(FIELD, x, y);
        let sum = 0;

        /// 4 point nbh
        if (ising_N === 1) {

            sum =   atGrid(FIELD, x + 1, y) + atGrid(FIELD, x - 1, y)
                +   atGrid(FIELD, x, y + 1) + atGrid(FIELD, x, y - 1);

        } 
        /// 8 point nbh
        else if (ising_N === 2) {

            sum =   atGrid(FIELD, x - 1, y - 1) +  atGrid(FIELD, x,     y - 1)
                +   atGrid(FIELD, x + 1, y - 1) +  atGrid(FIELD, x - 1, y    )
                +   atGrid(FIELD, x + 1, y    ) +  atGrid(FIELD, x - 1, y + 1)
                +   atGrid(FIELD, x,     y + 1) +  atGrid(FIELD, x + 1, y + 1); 
        } 
        /// other nbh
        else {
            for (let dx = -ising_N; dx <= ising_N; dx++) {
                for (let dy = -ising_N; dy <= ising_N; dy++) {
                    if (dx === 0 && dy === 0) continue;
                    sum += atGrid(FIELD, x + dx, y + dy) / Math.pow(Math.sqrt(dx * dx + dy * dy), 2);
                }
            }
        }

        const deltaE = 2 * spin * (ising_J * sum + ising_B);

        if ((deltaE < 0) || (Math.random() < Math.exp(-deltaE / ising_T))) 
        {
            FIELD[x][y] *= -1;
        }
    }
}


function sim_Run() {
    if (simModel === simCahn) {
        sim_cahnHilliardStep();
    } 
    else if (simModel === simIsing) {
        sim_isingStep();
    }

    sim_drawGrid();

    if (running) {
        simRequest = requestAnimationFrame(sim_Run);
    }
}

function sim_Start() {
    if (!running) {
        running = true;
        simRequest = requestAnimationFrame(sim_Run);
        document.getElementById('startBtn').disabled = true;
        document.getElementById('stopBtn').disabled  = false;
        document.getElementById('resetBtn').disabled = true;
    }
}

function sim_Stop() {
    if (running) {
        cancelAnimationFrame(simRequest);
        running = false;
        document.getElementById('startBtn').disabled = false;
        document.getElementById('stopBtn').disabled  = true;
        document.getElementById('resetBtn').disabled = false;
    }
}

function sim_Reset() {
    sim_Stop();
    sim_initModelGrid();
    sim_drawGrid();
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled  = true;
    document.getElementById('resetBtn').disabled = false;
}

function sim_UpdateGridSize(event) {
    sim_Stop();
    gridSize = parseInt(event.target.value);
    document.getElementById('gridSizeLabel').innerText = `Size: ${gridSize} x ${gridSize}`;
    sim_initModelGrid();
    sim_drawGrid();
    document.getElementById('startBtn').disabled = false;
    document.getElementById('stopBtn').disabled  = true;
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
    });

    paramDiv.appendChild(label);
    paramDiv.appendChild(slider);
}

function sim_SetInputParmeters() {
    const paramDiv = document.getElementById('parameters');
    paramDiv.innerHTML = '';

    if (simModel === simCahn) {
        
        sim_createSliderControl(paramDiv, {
            id: 'dxLabel',
            labelPrefix: 'Grid Spacing (Δx)',
            min: 0.1,
            max: 5,
            step: 0.01,
            value: cahn_dx,
            toFixed: 2,
            onChange: (val) => { cahn_dx = val; }
        });

        sim_createSliderControl(paramDiv, {
            id: 'dt',
            labelPrefix: 'Time sim_Run (Δt)',
            min: 0.001,
            max: 1,
            step: 0.001,
            value: cahn_dt,
            toFixed: 4,
            onChange: (val) => { cahn_dt = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'M',
            labelPrefix: 'Mobility (M)',
            min: 0.01,
            max: 1.0,
            step: 0.005,
            value: cahn_M,
            toFixed: 2,
            onChange: (val) => { cahn_M = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'k',
            labelPrefix: 'Interfacial Energy (κ)',
            min: 0.01,
            max: 1.0,
            step: 0.005,
            value: cahn_K,
            toFixed: 2,
            onChange: (val) => { cahn_K = val; }
        });
        
        sim_createSliderControl(paramDiv, {
            id: 'W',
            labelPrefix: 'Free Energy Coeff (A)',
            min: 0.1,
            max: 2.0,
            step: 0.1,
            value: cahn_W,
            toFixed: 2,
            onChange: (val) => { cahn_W = val; }
        });

    } else if (simModel === simIsing) {

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

        sim_createSliderControl(paramDiv, {
            id: 'N',
            labelPrefix: 'Neighborhood (N)',
            min: 1,
            max: 5,
            step: 1,
            value: ising_N,
            toFixed: 2,
            onChange: (val) => { ising_N = val; }
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

    dropdownContent.querySelectorAll('.dropdown-item').forEach(item => {
        item.addEventListener('click', function(event) {
            event.preventDefault();
            const selectedModel = this.getAttribute('data-model');

            simModel = selectedModel;

            menuToggle.firstChild.textContent = this.textContent + ' ';
            
            sim_Stop();
            sim_initModelGrid();
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
    document.getElementById('gridSizeInput').addEventListener('input', sim_UpdateGridSize);
    window.addEventListener('resize', sim_resizeCanvas);

    sim_resizeCanvas();
    sim_Stop();
};
