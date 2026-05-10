# ZigBio

ZigBio is a small static bio page that renders an interactive fluid simulation in Zig compiled to WebAssembly. The app is intentionally simple: one HTML file, one Zig simulation module, and no JavaScript framework.

## Current Decision

Keep this project as a static Zig/WASM page. The best upgrade path is not a framework rewrite; it is to strengthen the simulation logic, add repeatable build/test gates, and deploy the generated `zig-out/` directory to static hosting.

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

Deploy the contents of `zig-out/` to a static host such as GitHub Pages, Cloudflare Pages, Netlify, or Vercel. Configure `.wasm` files with `Content-Type: application/wasm` when the host supports custom MIME types. The page also includes a non-streaming fallback for hosts that serve WASM as `application/octet-stream`.
