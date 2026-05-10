# Project TODO

## Decision

Proceed with a static Zig/WASM deployment. The project is too small to justify a frontend framework or backend service. The near-term work should focus on simulation correctness, build reproducibility, polish, and deploy hygiene.

## P0 - Before Public Deployment

- [ ] Confirm the final public profile links in `index.html` before launch.
- [ ] Choose the deployment target: GitHub Pages, Cloudflare Pages, Netlify, or Vercel.
- [ ] Build with `zig build -Doptimize=ReleaseSmall` and deploy only the generated `zig-out/` contents.
- [ ] Verify the deployed host serves `fluid.wasm` correctly and the canvas works on mobile touch input.
- [ ] Add a short custom bio/name layer so the page is a real profile, not only a simulation demo.

## P1 - Zig Logic And Simulation Strengthening

- [ ] Make grid size, diffusion, viscosity, fade, and input force configurable at compile time.
- [ ] Add tests for boundary handling, velocity projection, density clamping, and one-step stability.
- [ ] Split the simulation into a reusable core module when the file grows beyond the current single-page scope.
- [ ] Add a lightweight benchmark step for `step()` so performance changes are measurable.

## P2 - Product Polish

- [ ] Add reduced-motion handling for users who prefer less animation.
- [ ] Tune mobile layout so profile links and controls never block touch interaction.
- [ ] Add a favicon, Open Graph metadata, and a preview image.
- [ ] Add a minimal CI workflow that runs `zig build test` and a release build.
