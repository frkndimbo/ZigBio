# TODO.md - ZigBio Floating Water World (Codex Edition, GBA-Authentic)

## System Prompt For Codex

You are building a GBA-authentic interactive portfolio page for a web browser.
The game logic, rendering, input state, camera, collision, navigation map, water
effect, and physics run entirely in Zig compiled to WebAssembly.

JavaScript is only allowed to:

- Load and instantiate the WASM module.
- Forward keyboard/touch input events into WASM input flags or axes.
- Read the final framebuffer from WASM and blit it to a `<canvas>`.
- Handle external navigation only after WASM exposes an active portal URL/state.

CSS is only for page chrome: centering the canvas, integer scaling, loading
screen, and optional non-game wrapper text. Zero game logic, sprite placement,
water animation, minimap drawing, or camera math in JS.

Read every section of this file before writing a single line of code.

## Hard Revision Decision

Performance and authenticity are the priority. Zig is the core game engine, not
just a helper for movement or fluid math. The current mixed JS renderer is an
intermediate prototype and must be migrated toward a Zig-owned framebuffer.

Visual direction: use the assets from `../Asset` and build a dense GBA-style
neighborhood inspired by Medabots (non-AX). Remove flood/swamp/floating themes.
Focus on streets, houses, fences, signs, trees, and clear route readability.

## The 90 Percent GBA Feel - Non-Negotiable Constraints

These constraints exist because GBA games look the way they do due to hardware
limits. Respecting them is what makes this feel authentic. Do not relax these
without explicit instruction.

### Screen

| Property | GBA Spec | This Project |
| --- | --- | --- |
| Resolution | 240 x 160 px | 240 x 160 internal render buffer |
| Display scaling | N/A | CSS `image-rendering: pixelated`; integer scale only, preferably 2x, 3x, or 4x |
| Color depth | 15-bit RGB555 | Scene palette max 256 colors |
| Refresh | About 60 Hz | `requestAnimationFrame`, target 60 fps |

The Zig WASM module maintains a `u16` framebuffer of `240 * 160` pixels in
RGB555 format. JS reads this buffer and blits it to canvas via `ImageData` at the
end of each frame.

Never render at native browser resolution and scale down. The pixelation must
come from the low internal resolution itself.

### Palette

- Maximum 256 colors in the scene palette, stored in WASM as `[256]u16`.
- All sprite and tile colors are palette-indexed: 4bpp or 8bpp style data.
- The water uses a dedicated 16-color sub-palette, matching GBA sprite palette
  bank behavior.
- Ripple effect is done by palette cycling or palette-index displacement, not by
  alpha blending.
- Ambient/background color is palette entry 0.

### Sprites

- Tile grid: 16 x 16 px per tile.
- Character sprite: 16 x 24 px, 3 tiles tall and 1 tile wide.
- Walk animation: 4 frames per direction, 4 directions, 16 frames total.
- Sprite sheet row-major order: down, up, left, right; 4 frames each.
- Houses: 32 x 32 px or 48 x 48 px sprites.
- Floating platforms, posts, signs, doors, roofs, and house details must reuse
  the assets from `../Asset` where possible.

### Audio

No audio in v1. Add only after all visual and performance tasks pass.

## Architecture - Mandatory, Do Not Deviate

### Zig Owns The Game

`fluid.zig` should evolve into a small GBA-like engine module, or be split into
focused Zig files if the build stays simple:

- `core.zig`: game loop entry points, fixed timestep, exported API.
- `framebuffer.zig`: RGB555 framebuffer, palette, clear/blit primitives.
- `tiles.zig`: tilemap renderer and palette-indexed tile data.
- `sprites.zig`: OAM-like sprite list, sprite drawing, animation frames.
- `water.zig`: existing fluid/ripple logic converted to palette-indexed water.
- `world.zig`: house positions, collision bounds, portal triggers, route data.
- `input.zig`: compact input state written by JS.
- `minimap.zig`: in-game map rendering so visitors do not get lost.

Keep the split pragmatic. Do not add module churn if one or two files are still
clearer, but the ownership boundaries above are mandatory.

### JavaScript Is A Cartridge Loader

`game.js` must become a thin loader:

- Instantiate `bin/fluid.wasm`.
- Create an `ImageData(240, 160)` scratch buffer.
- Convert each RGB555 pixel from WASM framebuffer to RGBA8888 only for blit.
- Forward input to exported Zig functions such as `setInput(...)`.
- Call exported `updateAndRender(dt_ms)` once per animation frame.
- Read exported portal state for navigation prompts.

Forbidden in JS:

- Drawing houses, player, water, particles, lights, paths, or minimap.
- Owning game positions or camera state.
- Procedural game-world rendering.
- Collision, proximity, pathfinding, animation timing, or map logic.

### CSS Is Page Chrome Only

CSS may style the browser shell, loading UI, integer scaling, and responsive
placement. CSS must not simulate game visuals or game HUD elements.

## Required Player Experience

- First screen shows the actual playable GBA-style world, not a landing page.
- Player starts near a useful central location with at least two reachable house
  directions visible or implied.
- Houses sit in readable neighborhood blocks connected by streets and signage.
- The map/minimap must show player position, house markers, and the current view
  rectangle or directional hints.
- Each house represents a profile portal: LinkedIn, GitHub, TikTok, Instagram,
  and Portfolio.
- Visitors should be able to understand where to go without reading long
  instructions.
- Mobile touch controls must work without covering the minimap or active portal
  prompt.

## Asset Rules

- Source assets live in `../Asset`.
- Build output may copy only the assets actually required for deployment.
- Prefer converting PNG art into palette-indexed Zig data during a build step or
  checked-in generated artifact.
- If a temporary runtime PNG loader remains during migration, document it as
  temporary and keep it out of the final architecture.
- Do not leave the scene plain. Use the house, roof, dock, sign, fence, item,
  and terrain tiles from the asset pack to build readable floating structures.

## Performance Budget

- Target: 60 fps on a mid-range mobile browser.
- Zig core update + render target: under 3 ms average per frame on desktop.
- JS blit/conversion target: under 2 ms average per frame on desktop.
- No per-frame heap allocation in Zig.
- No per-frame DOM changes except minimal status/navigation state.
- Fixed-size arrays for sprites, particles, ripples, houses, and map markers.
- Avoid alpha compositing as a core visual effect; prefer palette swaps, indexed
  tiles, and simple integer math.

## Implementation Plan

### Phase 1 - Lock The Contract

- [x] Update WASM exports around the final contract:
  - `init()`
  - `setInput(button_mask: u32, axis_x: f32, axis_y: f32)`
  - `updateAndRender(dt_ms: f32)`
  - `getFramebufferPtr() [*]u16`
  - `getPalettePtr() [*]u16`
  - `getActivePortal() i32`
  - `dismissActivePortal()`
- [x] Keep `getWorldWidth()` and `getWorldHeight()` only if useful for tests or
  debug; JS should not need them for rendering.
- [x] Add tests for framebuffer size, RGB555 color packing, input clamping, and
  active portal state.

### Phase 2 - Zig Framebuffer Renderer

- [x] Add RGB555 framebuffer storage in Zig: `[240 * 160]u16`.
- [x] Implement palette table with max 256 colors.
- [ ] Implement `clear`, `drawTile16`, `drawSprite`, and clipped blit routines.
- [x] Render water background in Zig using palette indices.
- [x] Render the player sprite in Zig.
- [x] Render floating houses and docks in Zig.
- [x] Render route markers and the minimap in Zig.

### Phase 3 - Asset Conversion

- [x] Websearch/reference mining >=100 aset GBA-style sebelum eksekusi (`docs/research/gba-asset-references.md`, 369 links).
- [x] Audit `../Asset/Tilemap/tilemap_packed.png` and selected `../Asset/Tiles`.
- [x] Select tiles for water props, floating bases, docks, roofs, doors, signs,
  and house bodies.
- [x] Convert selected tiles into palette-indexed data usable from Zig.
- [x] Keep generated data compact and deterministic.
- [x] Update `build.zig` so deployment installs only required runtime files.

### Phase 4 - Remove JS Renderer

- [x] Replace Canvas 2D world drawing in `game.js` with framebuffer blitting.
- [x] Delete JS house, player, water, route, minimap, and camera drawing logic.
- [x] Keep JS popup/navigation thin and state-driven from Zig.
- [x] Ensure no gameplay state lives in JS.

### Phase 5 - Game Feel And Navigation

- [ ] Tune player acceleration and walking speed for 240 x 160 readability.
- [ ] Add 4-direction walk animation with 4 frames per direction.
- [x] Remove flood/swamp visual dependency from the main world loop.
- [ ] Add clear minimap markers so visitors do not get lost.
- [ ] Tune house placement so the first viewport feels alive but not cluttered.
- [x] Apply per-portal logo pinpoints (LinkedIn, GitHub, TikTok, Instagram, Portfolio).
- [ ] Verify touch controls on mobile viewport.

### Phase 6 - Deployment Readiness

- [x] `rtk env ZIG_GLOBAL_CACHE_DIR=/tmp/bio-zig-global-cache ZIG_LOCAL_CACHE_DIR=/tmp/bio-zig-local-cache zig build test --summary all`
- [x] `rtk env ZIG_GLOBAL_CACHE_DIR=/tmp/bio-zig-global-cache ZIG_LOCAL_CACHE_DIR=/tmp/bio-zig-local-cache zig build -Doptimize=ReleaseSmall --summary all`
- [x] Browser screenshot validation, desktop and mobile.
- [x] Frame-time benchmark for Zig update/render and JS blit.
- [x] Add `wrangler.toml` for Cloudflare Pages with `pages_build_output_dir = "zig-out"`.
- [x] Document Cloudflare Pages deployment flow in `README.md`.
- [x] Execute first live deploy with authenticated Cloudflare account.
- [ ] Run `rtk graphify update .` after code edits. If graphify exits with the
  known `_os` NameError after writing output, note it instead of hiding it.
- [ ] Commit and push the finished revision to GitHub when validation passes.

## Current Prototype Gaps To Fix

- [x] JS still renders most of the game world. This violates the final
  architecture and must be removed.
- [x] Zig currently owns movement/proximity/fluid data, but not the final pixel
  framebuffer.
- [x] Water uses density/alpha-style overlay logic rather than true palette
  cycling.
- [x] Asset usage is visual but not yet converted into GBA-style indexed data.
- [x] The minimap exists conceptually but must be rendered by Zig in the final
  framebuffer.
- [x] The loading/page shell is acceptable, but game HUD must live inside the
  framebuffer.

## Definition Of Done

- [x] Browser page boots a WASM game at 240 x 160 internal resolution.
- [x] Zig owns game update, rendering, camera, collision, minimap, and
  portal state.
- [x] JS is limited to WASM loading, input forwarding, framebuffer blitting, and
  external navigation.
- [ ] Scene uses `../Asset` art for Medabots-inspired neighborhood housing.
- [ ] Visitors can navigate by looking at the in-game map.
- [x] ReleaseSmall build succeeds.
- [x] Unit tests pass.
- [x] Desktop and mobile visual checks pass.
- [ ] Graphify is updated after code edits.
- [ ] Changes are committed and pushed.
