(function () {
    const SCREEN_W = 240;
    const SCREEN_H = 160;
    const DEFAULT_WORLD_W = 900;
    const DEFAULT_WORLD_H = 620;
    const TILE = 16;
    const TILEMAP_SRC = 'Asset/Tilemap/tilemap_packed.png';

    const routeOrder = [0, 1, 2, 4, 3, 0];
    const houseLayouts = [
        { tiles: [[52, 53, 54, 55], [64, 65, 66, 67], [72, 73, 74, 75], [84, 85, 86, 87]], tileSize: 12, baseW: 60 },
        { tiles: [[48, 49, 50, 51], [60, 61, 62, 63], [76, 77, 78, 79], [88, 89, 90, 91]], tileSize: 12, baseW: 60 },
        { tiles: [[52, 53, 54], [64, 67, 68], [72, 74, 75], [84, 85, 87]], tileSize: 13, baseW: 56 },
        { tiles: [[48, 49, 50], [60, 63, 62], [76, 77, 79], [88, 89, 91]], tileSize: 13, baseW: 56 },
        { tiles: [[96, 97, 98, 99, 100], [108, 109, 110, 111, 112], [120, 121, 122, 123, 124]], tileSize: 12, baseW: 76 },
    ];

    const platforms = [
        { label: 'LinkedIn', short: 'IN', url: 'https://linkedin.com/in/frkndimbo', color: '#4aa3df' },
        { label: 'GitHub', short: 'GH', url: 'https://github.com/frkndimbo', color: '#dedede' },
        { label: 'TikTok', short: 'TK', url: 'https://www.tiktok.com/@frkndimbo', color: '#f05c74' },
        { label: 'Instagram', short: 'IG', url: 'https://www.instagram.com/frkndimbo', color: '#d980ff' },
        { label: 'Portfolio', short: 'PF', url: 'https://github.com/frkndimbo/ZigBio', color: '#f5c15b' },
    ];

    const prefersReducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const canvas = document.getElementById('game-canvas');
    const ctx = canvas.getContext('2d', { alpha: false });
    const loadingScreen = document.getElementById('loading-screen');
    const loadProgress = document.getElementById('load-progress');
    const popup = document.getElementById('game-popup');
    const popupIcon = document.getElementById('popup-icon');
    const popupTitle = document.getElementById('popup-title');
    const popupCopy = document.getElementById('popup-copy');
    const visitBtn = document.getElementById('visit-btn');
    const cancelBtn = document.getElementById('cancel-btn');
    const statusLabel = document.getElementById('status-label');
    const areaLabel = document.getElementById('area-label');
    const hintBanner = document.getElementById('hint-banner');
    const mapCanvas = document.getElementById('mini-map');
    const mapCtx = mapCanvas.getContext('2d', { alpha: true });
    const joystick = document.getElementById('joystick');
    const joystickKnob = document.getElementById('joystick-knob');
    const actionButton = document.getElementById('action-button');

    const state = {
        wasm: null,
        memory: null,
        worldW: DEFAULT_WORLD_W,
        worldH: DEFAULT_WORLD_H,
        gridSize: 0,
        houseCount: 0,
        rippleCount: 0,
        particleCount: 0,
        views: null,
        input: { x: 0, y: 0 },
        keys: new Set(),
        lastTime: 0,
        elapsed: 0,
        frame: 0,
        fps: 60,
        popupIndex: -1,
        tilesheet: null,
        houseSprites: [],
        backgroundCanvas: document.createElement('canvas'),
        waterCanvas: document.createElement('canvas'),
        waterCtx: null,
        waterData: null,
    };

    function clamp(value, min, max) {
        return Math.max(min, Math.min(max, value));
    }

    function setProgress(percent) {
        loadProgress.style.width = `${clamp(percent, 0, 100)}%`;
    }

    function loadImage(src) {
        return new Promise((resolve, reject) => {
            const image = new Image();
            image.onload = () => resolve(image);
            image.onerror = () => reject(new Error(`Unable to load ${src}`));
            image.src = src;
        });
    }

    async function loadTilesheet() {
        state.tilesheet = await loadImage(TILEMAP_SRC);
    }

    async function instantiateCore() {
        const imports = { env: { memory: new WebAssembly.Memory({ initial: 20 }) } };
        const candidates = ['bin/fluid.wasm', 'zig-out/bin/fluid.wasm'];
        let lastError = null;

        for (const url of candidates) {
            try {
                const response = await fetch(url);
                if (!response.ok) throw new Error(`${url} returned ${response.status}`);
                const contentType = response.headers.get('content-type') || '';
                const result = WebAssembly.instantiateStreaming && contentType.includes('application/wasm')
                    ? await WebAssembly.instantiateStreaming(Promise.resolve(response), imports)
                    : await WebAssembly.instantiate(await response.arrayBuffer(), imports);

                state.wasm = result.instance.exports;
                state.memory = state.wasm.memory || imports.env.memory;
                state.wasm.init();
                bindViews();
                return;
            } catch (error) {
                lastError = error;
            }
        }

        throw lastError || new Error('Unable to load Zig core');
    }

    function bindViews() {
        state.worldW = state.wasm.getWorldWidth ? state.wasm.getWorldWidth() : DEFAULT_WORLD_W;
        state.worldH = state.wasm.getWorldHeight ? state.wasm.getWorldHeight() : DEFAULT_WORLD_H;
        state.gridSize = state.wasm.getGridSize();
        state.houseCount = state.wasm.getHouseCount();
        state.rippleCount = state.wasm.getRippleCount();
        state.particleCount = state.wasm.getParticleCount();
        state.views = {
            houses: new Float32Array(state.memory.buffer, state.wasm.getHousePtr(), state.houseCount * 3),
            player: new Float32Array(state.memory.buffer, state.wasm.getPlayerPtr(), 6),
            camera: new Float32Array(state.memory.buffer, state.wasm.getCameraPtr(), 2),
            ripples: new Float32Array(state.memory.buffer, state.wasm.getRipplePtr(), state.rippleCount * 5),
            particles: new Float32Array(state.memory.buffer, state.wasm.getParticlePtr(), state.particleCount * 5),
            density: new Float32Array(state.memory.buffer, state.wasm.getDensityPtr(), (state.gridSize + 2) * (state.gridSize + 2)),
        };
    }

    function drawTileTo(target, tileIndex, dx, dy, size = TILE) {
        if (!state.tilesheet) return false;
        const sx = (tileIndex % 12) * TILE;
        const sy = Math.floor(tileIndex / 12) * TILE;
        target.imageSmoothingEnabled = false;
        target.drawImage(state.tilesheet, sx, sy, TILE, TILE, dx, dy, size, size);
        return true;
    }

    function drawRouteMarkers(target, from, to) {
        const dx = to.x - from.x;
        const dy = to.y - from.y;
        const distance = Math.hypot(dx, dy);
        const steps = Math.max(1, Math.floor(distance / 86));
        const markerTiles = [105, 107, 130, 131];

        for (let i = 1; i <= steps; i += 1) {
            const t = i / (steps + 1);
            const wave = Math.sin((from.x + to.y + i) * 0.08) * 8;
            const x = from.x + dx * t + wave;
            const y = from.y + dy * t - wave * 0.45;
            drawTileTo(target, markerTiles[i % markerTiles.length], x - 6, y - 6, 12);
        }
    }

    function buildBackground() {
        const bg = state.backgroundCanvas;
        bg.width = state.worldW;
        bg.height = state.worldH;
        const bctx = bg.getContext('2d', { alpha: false });
        bctx.fillStyle = '#277aa6';
        bctx.fillRect(0, 0, state.worldW, state.worldH);

        for (let y = 0; y < state.worldH; y += 8) {
            const shade = 14 + ((y / 8) % 4) * 4;
            bctx.fillStyle = `rgba(${25 + shade}, ${103 + shade}, ${142 + shade}, 0.22)`;
            bctx.fillRect(0, y, state.worldW, 3);
        }

        for (let y = 0; y < state.worldH; y += 48) {
            for (let x = (y / 3) % 64; x < state.worldW; x += 96) {
                bctx.fillStyle = 'rgba(194, 235, 222, 0.16)';
                bctx.fillRect(x, y + 18, 22, 2);
                bctx.fillRect(x + 7, y + 22, 16, 1);
            }
        }

        bctx.save();
        bctx.strokeStyle = 'rgba(246, 230, 183, 0.28)';
        bctx.lineWidth = 2;
        bctx.setLineDash([6, 10]);
        bctx.lineCap = 'square';
        bctx.beginPath();
        routeOrder.forEach((houseIndex, i) => {
            const house = houseAt(houseIndex);
            if (i === 0) bctx.moveTo(house.x, house.y + 13);
            else bctx.lineTo(house.x, house.y + 13);
        });
        bctx.stroke();
        bctx.restore();

        for (let i = 1; i < routeOrder.length; i += 1) {
            drawRouteMarkers(bctx, houseAt(routeOrder[i - 1]), houseAt(routeOrder[i]));
        }

        const edgeProps = [
            { tile: 93, x: 116, y: 118 },
            { tile: 95, x: state.worldW - 148, y: 118 },
            { tile: 105, x: 102, y: state.worldH - 134 },
            { tile: 107, x: state.worldW - 118, y: state.worldH - 92 },
        ];
        edgeProps.forEach((prop) => drawTileTo(bctx, prop.tile, prop.x, prop.y, 14));
    }

    function buildHouseSprites() {
        state.houseSprites = houseLayouts.map((layout) => {
            const columns = layout.tiles.reduce((max, row) => Math.max(max, row.length), 0);
            const canvas = document.createElement('canvas');
            canvas.width = columns * layout.tileSize;
            canvas.height = layout.tiles.length * layout.tileSize;
            const hctx = canvas.getContext('2d', { alpha: true });
            layout.tiles.forEach((row, rowIndex) => {
                row.forEach((tileIndex, columnIndex) => {
                    drawTileTo(hctx, tileIndex, columnIndex * layout.tileSize, rowIndex * layout.tileSize, layout.tileSize);
                });
            });
            return {
                canvas,
                width: canvas.width,
                height: canvas.height,
                baseW: layout.baseW,
                anchorY: canvas.height - 8,
                markerOffset: canvas.height + 10,
            };
        });
    }

    function prepareWaterCanvas() {
        state.waterCanvas.width = state.gridSize;
        state.waterCanvas.height = state.gridSize;
        state.waterCtx = state.waterCanvas.getContext('2d', { alpha: true });
        state.waterData = state.waterCtx.createImageData(state.gridSize, state.gridSize);
    }

    function inputVector() {
        let x = state.input.x;
        let y = state.input.y;
        if (state.keys.has('ArrowLeft') || state.keys.has('KeyA')) x -= 1;
        if (state.keys.has('ArrowRight') || state.keys.has('KeyD')) x += 1;
        if (state.keys.has('ArrowUp') || state.keys.has('KeyW')) y -= 1;
        if (state.keys.has('ArrowDown') || state.keys.has('KeyS')) y += 1;

        const length = Math.hypot(x, y);
        if (length > 1) {
            x /= length;
            y /= length;
        }
        return { x, y };
    }

    function updatePopup() {
        const nextIndex = state.wasm.getActiveHouse();
        if (nextIndex === state.popupIndex) return;
        state.popupIndex = nextIndex;

        if (nextIndex < 0) {
            popup.hidden = true;
            areaLabel.textContent = 'FLOAT MAP';
            return;
        }

        const platform = platforms[nextIndex];
        popupIcon.textContent = platform.short;
        popupIcon.style.backgroundColor = platform.color;
        popupTitle.textContent = platform.label;
        popupCopy.textContent = 'Mau berkunjung?';
        areaLabel.textContent = `${platform.label} HOUSE`;
        popup.hidden = false;
        window.ZigBioPreloader.prefetchUrl(platform.url);
    }

    function update(dt) {
        const input = inputVector();
        state.wasm.gameUpdate(dt, input.x, input.y, prefersReducedMotion ? 1 : 0);
        updatePopup();
    }

    function updateWaterImage() {
        if (state.frame % 2 !== 0) return;
        const pixels = state.waterData.data;
        const density = state.views.density;
        const size = state.gridSize;

        for (let y = 1; y <= size; y += 1) {
            for (let x = 1; x <= size; x += 1) {
                const value = clamp(density[x + (size + 2) * y] * 0.54, 0, 180);
                const p = ((y - 1) * size + (x - 1)) * 4;
                pixels[p] = 115;
                pixels[p + 1] = 199;
                pixels[p + 2] = 188;
                pixels[p + 3] = value;
            }
        }
        state.waterCtx.putImageData(state.waterData, 0, 0);
    }

    function drawBackground(cameraX, cameraY) {
        ctx.drawImage(state.backgroundCanvas, cameraX, cameraY, SCREEN_W, SCREEN_H, 0, 0, SCREEN_W, SCREEN_H);
    }

    function drawWater(cameraX, cameraY) {
        updateWaterImage();
        const scaleX = state.gridSize / state.worldW;
        const scaleY = state.gridSize / state.worldH;
        ctx.save();
        ctx.globalAlpha = 0.68;
        ctx.drawImage(
            state.waterCanvas,
            cameraX * scaleX,
            cameraY * scaleY,
            SCREEN_W * scaleX,
            SCREEN_H * scaleY,
            0,
            0,
            SCREEN_W,
            SCREEN_H,
        );
        ctx.restore();
    }

    function houseAt(index) {
        const base = index * 3;
        return {
            x: state.views.houses[base],
            y: state.views.houses[base + 1],
            radius: state.views.houses[base + 2],
            platform: platforms[index],
        };
    }

    function drawLights(cameraX, cameraY) {
        const px = state.views.player[0];
        const py = state.views.player[1];
        ctx.save();
        ctx.globalCompositeOperation = 'screen';
        for (let i = 0; i < state.houseCount; i += 1) {
            const house = houseAt(i);
            const sx = house.x - cameraX;
            const sy = house.y - cameraY - 4;
            const distance = Math.hypot(px - house.x, py - house.y);
            const proximity = clamp(1 - distance / 170, 0, 1);
            const flicker = prefersReducedMotion ? 0 : Math.sin(state.elapsed * 5.2 + house.x) * 0.04;
            const radius = 42 + proximity * 32;
            const gradient = ctx.createRadialGradient(sx, sy, 3, sx, sy, radius);
            gradient.addColorStop(0, `rgba(245, 193, 91, ${0.46 + proximity * 0.34 + flicker})`);
            gradient.addColorStop(1, 'rgba(245, 193, 91, 0)');
            ctx.fillStyle = gradient;
            ctx.fillRect(sx - radius, sy - radius, radius * 2, radius * 2);
        }
        ctx.restore();
    }

    function drawRipples(cameraX, cameraY) {
        if (prefersReducedMotion) return;
        const ripples = state.views.ripples;
        ctx.save();
        ctx.strokeStyle = '#b8f3d4';
        ctx.lineWidth = 1;
        for (let i = 0; i < state.rippleCount; i += 1) {
            const base = i * 5;
            if (ripples[base + 4] < 0.5) continue;
            ctx.globalAlpha = ripples[base + 3];
            ctx.beginPath();
            ctx.ellipse(ripples[base] - cameraX, ripples[base + 1] - cameraY, ripples[base + 2] * 1.3, ripples[base + 2] * 0.42, 0, 0, Math.PI * 2);
            ctx.stroke();
        }
        ctx.restore();

        const particles = state.views.particles;
        ctx.save();
        ctx.fillStyle = '#d7ffdf';
        for (let i = 0; i < state.particleCount; i += 1) {
            const base = i * 5;
            const life = particles[base + 4];
            if (life <= 0) continue;
            ctx.globalAlpha = clamp(life * 2.7, 0, 1);
            ctx.fillRect(Math.round(particles[base] - cameraX), Math.round(particles[base + 1] - cameraY), 1, 1);
        }
        ctx.restore();
    }

    function drawFloatingBase(x, y, width) {
        const half = width / 2;
        ctx.save();
        ctx.fillStyle = 'rgba(8, 24, 36, 0.44)';
        ctx.beginPath();
        ctx.ellipse(x, y + 16, half + 7, 13, 0, 0, Math.PI * 2);
        ctx.fill();

        ctx.fillStyle = '#7c4a31';
        ctx.fillRect(x - half, y + 4, width, 12);
        ctx.fillStyle = '#b77a48';
        for (let offset = -half + 4; offset < half - 4; offset += 12) {
            ctx.fillRect(Math.round(x + offset), y + 5, 8, 9);
        }
        ctx.fillStyle = '#4a2a27';
        ctx.fillRect(x - half - 2, y + 1, width + 4, 4);
        ctx.fillRect(x - half - 2, y + 16, width + 4, 4);
        ctx.restore();
    }

    function drawHouse(house, index, cameraX, cameraY) {
        const x = Math.round(house.x - cameraX);
        const y = Math.round(house.y - cameraY);

        const bob = prefersReducedMotion ? 0 : Math.sin(state.elapsed * 1.8 + index) * 1.5;
        const sprite = state.houseSprites.length > 0 ? state.houseSprites[index % state.houseSprites.length] : null;
        drawFloatingBase(x, y + 8 + bob, sprite ? sprite.baseW : 60);

        if (sprite) {
            ctx.drawImage(sprite.canvas, Math.round(x - sprite.width / 2), Math.round(y - sprite.anchorY + bob));
            return;
        }

        ctx.fillStyle = '#11191e';
        ctx.fillRect(x - 24, y - 8 + bob, 48, 18);
        ctx.fillStyle = '#7c4538';
        ctx.fillRect(x - 18, y - 24 + bob, 36, 26);
        ctx.fillStyle = '#3b2430';
        ctx.beginPath();
        ctx.moveTo(x - 23, y - 22 + bob);
        ctx.lineTo(x, y - 42 + bob);
        ctx.lineTo(x + 23, y - 22 + bob);
        ctx.closePath();
        ctx.fill();
    }

    function drawPinpoint(house, index, cameraX, cameraY) {
        const sprite = state.houseSprites.length > 0 ? state.houseSprites[index % state.houseSprites.length] : null;
        const markerOffset = sprite ? sprite.markerOffset : 46;
        const bob = prefersReducedMotion ? 0 : Math.sin(state.elapsed * 4 + house.x * 0.01) * 3;
        const x = Math.round(house.x - cameraX);
        const y = Math.round(house.y - cameraY - markerOffset + bob);

        ctx.fillStyle = '#101820';
        ctx.fillRect(x - 13, y - 13, 26, 20);
        ctx.fillStyle = house.platform.color;
        ctx.fillRect(x - 10, y - 10, 20, 14);
        ctx.fillStyle = '#101820';
        ctx.beginPath();
        ctx.moveTo(x - 4, y + 6);
        ctx.lineTo(x + 4, y + 6);
        ctx.lineTo(x, y + 13);
        ctx.closePath();
        ctx.fill();
        ctx.fillStyle = '#101820';
        ctx.font = 'bold 8px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(house.platform.short, x, y - 3);
    }

    function drawPlayer(cameraX, cameraY) {
        const player = state.views.player;
        const x = Math.round(player[0] - cameraX);
        const y = Math.round(player[1] - cameraY);
        const dir = player[4];
        const moving = player[5] > 0.5;
        const bob = moving && !prefersReducedMotion ? Math.sin(state.elapsed * 16) * 2 : 0;

        ctx.fillStyle = 'rgba(0, 0, 0, 0.28)';
        ctx.fillRect(x - 8, y + 8, 16, 4);
        ctx.fillStyle = '#242f37';
        ctx.fillRect(x - 5, y - 14 + bob, 10, 7);
        ctx.fillStyle = '#f0c08a';
        ctx.fillRect(x - 4, y - 21 + bob, 8, 7);
        ctx.fillStyle = '#263b42';
        ctx.fillRect(x - 6, y - 25 + bob, 12, 5);
        ctx.fillStyle = '#f5c15b';
        ctx.fillRect(x - 6, y - 7 + bob, 5, 11);
        ctx.fillRect(x + 1, y - 7 - bob, 5, 11);
        ctx.fillStyle = '#101820';
        ctx.fillRect(x - 6, y + 4 + bob, 5, 3);
        ctx.fillRect(x + 1, y + 4 - bob, 5, 3);

        if (dir === 2) ctx.fillRect(x - 5, y - 18 + bob, 2, 2);
        else if (dir === 3) ctx.fillRect(x + 3, y - 18 + bob, 2, 2);
        else {
            ctx.fillRect(x - 3, y - 18 + bob, 2, 2);
            ctx.fillRect(x + 2, y - 18 + bob, 2, 2);
        }
    }

    function drawVignette() {
        const gradient = ctx.createRadialGradient(SCREEN_W / 2, SCREEN_H / 2, 30, SCREEN_W / 2, SCREEN_H / 2, 150);
        gradient.addColorStop(0, 'rgba(0, 0, 0, 0)');
        gradient.addColorStop(1, 'rgba(0, 0, 0, 0.28)');
        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, SCREEN_W, SCREEN_H);
    }

    function drawMiniMap(cameraX, cameraY) {
        const scaleX = mapCanvas.width / state.worldW;
        const scaleY = mapCanvas.height / state.worldH;
        mapCtx.clearRect(0, 0, mapCanvas.width, mapCanvas.height);
        mapCtx.fillStyle = '#14384e';
        mapCtx.fillRect(0, 0, mapCanvas.width, mapCanvas.height);

        mapCtx.fillStyle = 'rgba(115, 199, 188, 0.22)';
        for (let y = 5; y < mapCanvas.height; y += 9) {
            mapCtx.fillRect(4, y, mapCanvas.width - 8, 1);
        }

        for (let i = 0; i < state.houseCount; i += 1) {
            const house = houseAt(i);
            mapCtx.fillStyle = house.platform.color;
            mapCtx.fillRect(Math.round(house.x * scaleX) - 2, Math.round(house.y * scaleY) - 2, 4, 4);
        }

        mapCtx.strokeStyle = '#f6e6b7';
        mapCtx.lineWidth = 1;
        mapCtx.strokeRect(cameraX * scaleX, cameraY * scaleY, SCREEN_W * scaleX, SCREEN_H * scaleY);

        mapCtx.fillStyle = '#fff4cf';
        mapCtx.fillRect(Math.round(state.views.player[0] * scaleX) - 1, Math.round(state.views.player[1] * scaleY) - 1, 3, 3);
    }

    function render() {
        const cameraX = Math.round(state.views.camera[0]);
        const cameraY = Math.round(state.views.camera[1]);
        ctx.imageSmoothingEnabled = false;
        drawBackground(cameraX, cameraY);
        drawWater(cameraX, cameraY);
        drawRipples(cameraX, cameraY);
        drawLights(cameraX, cameraY);

        const drawables = [{ type: 'player', y: state.views.player[1] }];
        for (let i = 0; i < state.houseCount; i += 1) {
            const house = houseAt(i);
            drawables.push({ type: 'house', house, index: i, y: house.y });
        }
        drawables.sort((a, b) => a.y - b.y);
        drawables.forEach((item) => {
            if (item.type === 'player') drawPlayer(cameraX, cameraY);
            else drawHouse(item.house, item.index, cameraX, cameraY);
        });

        for (let i = 0; i < state.houseCount; i += 1) drawPinpoint(houseAt(i), i, cameraX, cameraY);
        drawMiniMap(cameraX, cameraY);
        drawVignette();
    }

    function tick(now) {
        if (!state.lastTime) state.lastTime = now;
        const dt = Math.min((now - state.lastTime) / 1000, 0.05);
        state.lastTime = now;
        state.elapsed += dt;
        state.frame += 1;
        state.fps = state.fps * 0.92 + (1 / Math.max(dt, 0.001)) * 0.08;

        update(dt);
        render();
        if (state.frame % 30 === 0) statusLabel.textContent = `${Math.round(state.fps)} FPS`;
        requestAnimationFrame(tick);
    }

    function visitActiveHouse() {
        const index = state.wasm.getActiveHouse();
        if (index < 0) return;
        window.ZigBioPreloader.prefetchAndNavigate(platforms[index].url);
    }

    function hidePopup() {
        state.wasm.dismissActiveHouse();
        state.popupIndex = -1;
        popup.hidden = true;
        areaLabel.textContent = 'FLOAT MAP';
    }

    function setupKeyboard() {
        window.addEventListener('keydown', (event) => {
            if (['ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight', 'Space'].includes(event.code)) event.preventDefault();
            state.keys.add(event.code);
            if (event.code === 'Enter' || event.code === 'Space') visitActiveHouse();
            if (event.code === 'Escape') hidePopup();
        });
        window.addEventListener('keyup', (event) => state.keys.delete(event.code));
    }

    function setupPopup() {
        visitBtn.addEventListener('click', visitActiveHouse);
        cancelBtn.addEventListener('click', hidePopup);
    }

    function setupJoystick() {
        let activePointer = null;
        const resetJoystick = () => {
            activePointer = null;
            state.input.x = 0;
            state.input.y = 0;
            joystickKnob.style.transform = 'translate(-50%, -50%)';
        };

        joystick.addEventListener('pointerdown', (event) => {
            activePointer = event.pointerId;
            joystick.setPointerCapture(activePointer);
        });
        joystick.addEventListener('pointermove', (event) => {
            if (event.pointerId !== activePointer) return;
            const rect = joystick.getBoundingClientRect();
            const centerX = rect.left + rect.width / 2;
            const centerY = rect.top + rect.height / 2;
            const dx = event.clientX - centerX;
            const dy = event.clientY - centerY;
            const max = rect.width * 0.34;
            const length = Math.min(Math.hypot(dx, dy), max);
            const angle = Math.atan2(dy, dx);
            const x = Math.cos(angle) * length;
            const y = Math.sin(angle) * length;
            state.input.x = x / max;
            state.input.y = y / max;
            joystickKnob.style.transform = `translate(calc(-50% + ${x}px), calc(-50% + ${y}px))`;
        });
        joystick.addEventListener('pointerup', resetJoystick);
        joystick.addEventListener('pointercancel', resetJoystick);
        actionButton.addEventListener('click', visitActiveHouse);
    }

    async function boot() {
        canvas.width = SCREEN_W;
        canvas.height = SCREEN_H;
        setProgress(12);
        await loadTilesheet();
        setProgress(42);
        await instantiateCore();
        setProgress(74);
        buildHouseSprites();
        buildBackground();
        prepareWaterCanvas();
        setProgress(100);
        setupKeyboard();
        setupPopup();
        setupJoystick();
        statusLabel.textContent = '60 FPS';
        render();
        loadingScreen.classList.add('is-done');
        setTimeout(() => hintBanner.classList.add('is-hidden'), 3200);
        requestAnimationFrame(tick);
    }

    boot().catch(() => {
        statusLabel.textContent = 'BOOT ERR';
        popup.hidden = true;
        setProgress(100);
    });
})();
