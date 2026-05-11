# ZigBio

ZigBio is a static GBA-style bio game. The game core runs in Zig compiled to WebAssembly and now owns movement, camera, proximity, minimap, water rendering, and the RGB555 framebuffer. JavaScript is limited to input forwarding, framebuffer blitting, and navigation UI.

## Current Decision

Keep this project static and performance-first. Do not add a frontend framework or backend service. The browser should load a small WASM core and static files from `zig-out/`.

## Requirements

- Zig 0.16.0 or newer compatible 0.16.x build system
- A static web server for local preview, because browser `fetch()` cannot reliably load WASM from `file://`

## Build

```sh
zig build -Doptimize=ReleaseSmall
```

For this local Codex environment, use writable Zig caches:

```sh
rtk env ZIG_GLOBAL_CACHE_DIR=/tmp/bio-zig-global-cache ZIG_LOCAL_CACHE_DIR=/tmp/bio-zig-local-cache zig build -Doptimize=ReleaseSmall
```

The deployable output is written to:

```text
zig-out/
├── index.html
├── game.js
├── preloader.js
├── style.css
└── bin/fluid.wasm
```

## Test

```sh
zig build test
```

## Local Preview

After building, serve `zig-out/` with any static server:

```sh
python3 -m http.server 4173 --directory zig-out
```

Then open `http://127.0.0.1:4173/`.

## Deployment

Deploy the contents of `zig-out/` to a static host such as GitHub Pages, Cloudflare Pages, Netlify, or Vercel. Configure `.wasm` files with `Content-Type: application/wasm` when the host supports custom MIME types. The loader includes a non-streaming fallback for hosts that serve WASM as `application/octet-stream`.
