# Graph Report - fluid-sim  (2026-05-11)

## Corpus Check
- 6 files · ~11,042 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 131 nodes · 243 edges · 15 communities detected
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]

## God Nodes (most connected - your core abstractions)
1. `render()` - 11 edges
2. `boot()` - 11 edges
3. `render()` - 11 edges
4. `boot()` - 11 edges
5. `clamp()` - 8 edges
6. `idx()` - 8 edges
7. `gameUpdate()` - 8 edges
8. `updatePlayer()` - 7 edges
9. `syncViews()` - 6 edges
10. `project()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `setPlayerPosition()` --calls--> `clamp()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 3 → community 6_
- `advect()` --calls--> `clamp()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 3 → community 5_
- `drawLights()` --calls--> `clamp()`  [EXTRACTED]
  game.js → game.js  _Bridges community 10 → community 8_
- `boot()` --calls--> `setProgress()`  [EXTRACTED]
  game.js → game.js  _Bridges community 10 → community 7_
- `boot()` --calls--> `instantiateCore()`  [EXTRACTED]
  game.js → game.js  _Bridges community 2 → community 7_

## Communities

### Community 0 - "Community 0"
Cohesion: 0.09
Nodes (6): House, Particle, Player, Ripple, Screen, World

### Community 1 - "Community 1"
Cohesion: 0.23
Nodes (9): buildBackground(), drawFloatingBase(), drawHouse(), drawRouteMarkers(), drawTileTo(), inputVector(), tick(), update() (+1 more)

### Community 2 - "Community 2"
Cohesion: 0.24
Nodes (8): bindViews(), drawFloatingBase(), drawHouse(), inputVector(), instantiateCore(), tick(), update(), updatePopup()

### Community 3 - "Community 3"
Cohesion: 0.27
Nodes (11): addDensity(), clamp(), distance(), disturbFluid(), gameUpdate(), spawnRipple(), spawnSplash(), updateCamera() (+3 more)

### Community 4 - "Community 4"
Cohesion: 0.2
Nodes (10): bindViews(), boot(), buildHouseSprites(), instantiateCore(), loadImage(), loadTilesheet(), prepareWaterCanvas(), setupJoystick() (+2 more)

### Community 5 - "Community 5"
Cohesion: 0.42
Nodes (9): addVelocity(), advect(), diffuse(), divergenceAt(), idx(), linSolve(), project(), setBoundary() (+1 more)

### Community 6 - "Community 6"
Cohesion: 0.22
Nodes (9): init(), resetEffects(), resetFluid(), setPlayerPosition(), syncHouseView(), syncParticleView(), syncPlayerView(), syncRippleView() (+1 more)

### Community 7 - "Community 7"
Cohesion: 0.25
Nodes (8): boot(), buildHouseSprites(), loadImage(), loadTilesheet(), prepareWaterCanvas(), setupJoystick(), setupKeyboard(), setupPopup()

### Community 8 - "Community 8"
Cohesion: 0.38
Nodes (7): drawBackground(), drawLights(), drawMiniMap(), drawPinpoint(), drawVignette(), houseAt(), render()

### Community 9 - "Community 9"
Cohesion: 0.38
Nodes (7): drawBackground(), drawLights(), drawMiniMap(), drawPinpoint(), drawVignette(), houseAt(), render()

### Community 10 - "Community 10"
Cohesion: 0.4
Nodes (5): clamp(), drawRipples(), drawWater(), setProgress(), updateWaterImage()

### Community 11 - "Community 11"
Cohesion: 0.4
Nodes (5): clamp(), drawRipples(), drawWater(), setProgress(), updateWaterImage()

### Community 12 - "Community 12"
Cohesion: 0.83
Nodes (3): prefetchAndNavigate(), prefetchUrl(), waitForDestination()

### Community 13 - "Community 13"
Cohesion: 0.83
Nodes (3): prefetchAndNavigate(), prefetchUrl(), waitForDestination()

### Community 14 - "Community 14"
Cohesion: 0.67
Nodes (3): buildBackground(), drawRouteMarkers(), drawTileTo()

## Knowledge Gaps
- **6 isolated node(s):** `Screen`, `World`, `House`, `Player`, `Ripple` (+1 more)
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `boot()` connect `Community 7` to `Community 8`, `Community 2`, `Community 10`, `Community 14`?**
  _High betweenness centrality (0.003) - this node is a cross-community bridge._
- **Why does `boot()` connect `Community 4` to `Community 1`, `Community 11`, `Community 9`?**
  _High betweenness centrality (0.003) - this node is a cross-community bridge._
- **Why does `render()` connect `Community 8` to `Community 2`, `Community 10`, `Community 7`?**
  _High betweenness centrality (0.003) - this node is a cross-community bridge._
- **What connects `Screen`, `World`, `House` to the rest of the system?**
  _6 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.09 - nodes in this community are weakly interconnected._