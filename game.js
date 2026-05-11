(function () {
    const SCREEN_W = 240;
    const SCREEN_H = 160;

    const BUTTON_REDUCED_MOTION = 1 << 2;

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
    const joystick = document.getElementById('joystick');
    const joystickKnob = document.getElementById('joystick-knob');
    const actionButton = document.getElementById('action-button');

    const state = {
        wasm: null,
        memory: null,
        framebuffer: null,
        palette: null,
        rgbLookup: new Uint32Array(32768),
        imageData: null,
        rgba32: null,
        input: { x: 0, y: 0 },
        keys: new Set(),
        lastTime: 0,
        frame: 0,
        fps: 60,
        popupIndex: -1,
    };

    function clamp(value, min, max) {
        return Math.max(min, Math.min(max, value));
    }

    function setProgress(percent) {
        loadProgress.style.width = `${clamp(percent, 0, 100)}%`;
    }

    function buildRgbLookup() {
        for (let color = 0; color < state.rgbLookup.length; color += 1) {
            const r5 = color & 0x1f;
            const g5 = (color >>> 5) & 0x1f;
            const b5 = (color >>> 10) & 0x1f;
            const r8 = (r5 << 3) | (r5 >>> 2);
            const g8 = (g5 << 3) | (g5 >>> 2);
            const b8 = (b5 << 3) | (b5 >>> 2);
            state.rgbLookup[color] = (255 << 24) | (b8 << 16) | (g8 << 8) | r8;
        }
    }

    function activePortalIndex() {
        if (state.wasm.getActivePortal) return state.wasm.getActivePortal();
        return state.wasm.getActiveHouse();
    }

    function dismissPortal() {
        if (state.wasm.dismissActivePortal) state.wasm.dismissActivePortal();
        else state.wasm.dismissActiveHouse();
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
        const framebufferLen = state.wasm.getFramebufferLen ? state.wasm.getFramebufferLen() : SCREEN_W * SCREEN_H;
        state.framebuffer = new Uint16Array(state.memory.buffer, state.wasm.getFramebufferPtr(), framebufferLen);
        if (state.wasm.getPalettePtr) {
            state.palette = new Uint16Array(state.memory.buffer, state.wasm.getPalettePtr(), 256);
        }
    }

    function blitFrame() {
        const src = state.framebuffer;
        const dst = state.rgba32;
        const lut = state.rgbLookup;
        for (let i = 0; i < src.length; i += 1) {
            dst[i] = lut[src[i] & 0x7fff];
        }
        ctx.putImageData(state.imageData, 0, 0);
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
        const nextIndex = activePortalIndex();
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

    function stepCore(dtMs) {
        const input = inputVector();
        let buttonMask = 0;
        if (prefersReducedMotion) buttonMask |= BUTTON_REDUCED_MOTION;
        state.wasm.setInput(buttonMask, input.x, input.y);
        state.wasm.updateAndRender(dtMs);
    }

    function tick(now) {
        if (!state.lastTime) state.lastTime = now;
        const dtMs = Math.min(now - state.lastTime, 50);
        state.lastTime = now;
        state.frame += 1;
        state.fps = state.fps * 0.92 + (1000 / Math.max(dtMs, 1)) * 0.08;

        stepCore(dtMs);
        updatePopup();
        blitFrame();

        if (state.frame % 30 === 0) statusLabel.textContent = `${Math.round(state.fps)} FPS`;
        requestAnimationFrame(tick);
    }

    function visitActiveHouse() {
        const index = activePortalIndex();
        if (index < 0) return;
        window.ZigBioPreloader.prefetchAndNavigate(platforms[index].url);
    }

    function hidePopup() {
        dismissPortal();
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
        ctx.imageSmoothingEnabled = false;

        state.imageData = ctx.createImageData(SCREEN_W, SCREEN_H);
        state.rgba32 = new Uint32Array(state.imageData.data.buffer);
        buildRgbLookup();

        setProgress(22);
        await instantiateCore();
        setProgress(90);

        setupKeyboard();
        setupPopup();
        setupJoystick();

        statusLabel.textContent = '60 FPS';
        blitFrame();
        setProgress(100);

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
