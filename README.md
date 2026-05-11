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

Asset tiles for houses are converted to palette-indexed Zig data in
`asset_tiles.zig`. Regenerate them from `../Asset/Tiles` with:

```sh
rtk python3 tools/generate_asset_tiles.py
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

### Cloudflare Pages (Direct Upload)

1. Build release artifacts:

   ```sh
   rtk env ZIG_GLOBAL_CACHE_DIR=/tmp/bio-zig-global-cache ZIG_LOCAL_CACHE_DIR=/tmp/bio-zig-local-cache zig build -Doptimize=ReleaseSmall
   ```

2. Authenticate Wrangler once:

   ```sh
   rtk bunx wrangler login
   ```

3. Create a Pages project once:

   ```sh
   rtk bunx wrangler pages project create zigbio --production-branch main
   ```

4. Deploy the generated static output:

   ```sh
   rtk bunx wrangler pages deploy zig-out --project-name zigbio --branch main
   ```

`wrangler.toml` is included with `pages_build_output_dir = "zig-out"` for this
project. If your project name differs, update `name` in `wrangler.toml` and the
deploy command accordingly.

### Other Static Hosts

You can also deploy `zig-out/` to GitHub Pages, Netlify, or Vercel. The loader
already includes a fallback for hosts that return WASM as
`application/octet-stream`.
