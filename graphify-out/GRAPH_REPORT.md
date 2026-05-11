# Graph Report - fluid-sim  (2026-05-11)

## Corpus Check
- 8 files · ~47,846 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 134 nodes · 243 edges · 14 communities detected
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

## God Nodes (most connected - your core abstractions)
1. `renderFrame()` - 11 edges
2. `clamp()` - 10 edges
3. `fillRect()` - 9 edges
4. `idx()` - 8 edges
5. `init()` - 8 edges
6. `boot()` - 8 edges
7. `boot()` - 8 edges
8. `stepGame()` - 7 edges
9. `drawHouse()` - 6 edges
10. `syncViews()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `setPlayerPosition()` --calls--> `clamp()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 3 → community 8_
- `disturbFluid()` --calls--> `clamp()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 3 → community 4_
- `drawMapTile()` --calls--> `drawAssetTile16()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 9 → community 7_
- `drawHouse()` --calls--> `drawAssetTile16()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 9 → community 5_
- `drawTownMap()` --calls--> `drawMapTile()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 7 → community 5_

## Communities

### Community 0 - "Community 0"
Cohesion: 0.07
Nodes (7): House, HouseSprite, Particle, Player, Ripple, Screen, World

### Community 1 - "Community 1"
Cohesion: 0.2
Nodes (18): activePortalIndex(), bindViews(), blitFrame(), boot(), buildRgbLookup(), clamp(), dismissPortal(), hidePopup() (+10 more)

### Community 2 - "Community 2"
Cohesion: 0.2
Nodes (18): activePortalIndex(), bindViews(), blitFrame(), boot(), buildRgbLookup(), clamp(), dismissPortal(), hidePopup() (+10 more)

### Community 3 - "Community 3"
Cohesion: 0.25
Nodes (11): clamp(), distance(), gameUpdate(), reducedMotionEnabled(), setInput(), spawnRipple(), stepGame(), updateAndRender() (+3 more)

### Community 4 - "Community 4"
Cohesion: 0.33
Nodes (11): addDensity(), addVelocity(), advect(), diffuse(), disturbFluid(), divergenceAt(), idx(), linSolve() (+3 more)

### Community 5 - "Community 5"
Cohesion: 0.4
Nodes (10): clampI(), drawFrameRect(), drawHouse(), drawMiniMap(), drawPlayer(), drawPortalBadge(), drawRouteDots(), drawTownMap() (+2 more)

### Community 6 - "Community 6"
Cohesion: 0.33
Nodes (6): init(), initPalette(), resetEffects(), resetFluid(), rgb555(), syncHouseView()

### Community 7 - "Community 7"
Cohesion: 0.5
Nodes (5): drawMapTile(), fillMapRect(), initTownMap(), mapOffset(), setMapTile()

### Community 8 - "Community 8"
Cohesion: 0.4
Nodes (5): setPlayerPosition(), syncParticleView(), syncPlayerView(), syncRippleView(), syncViews()

### Community 9 - "Community 9"
Cohesion: 0.5
Nodes (4): drawAssetTile16(), fbOffset(), findAssetTile(), setPixel()

### Community 10 - "Community 10"
Cohesion: 0.83
Nodes (3): prefetchAndNavigate(), prefetchUrl(), waitForDestination()

### Community 11 - "Community 11"
Cohesion: 0.83
Nodes (3): prefetchAndNavigate(), prefetchUrl(), waitForDestination()

### Community 12 - "Community 12"
Cohesion: 1.0
Nodes (2): main(), rgb555()

### Community 13 - "Community 13"
Cohesion: 1.0
Nodes (2): dismissActiveHouse(), dismissActivePortal()

## Knowledge Gaps
- **7 isolated node(s):** `Screen`, `World`, `House`, `Player`, `Ripple` (+2 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 12`** (3 nodes): `main()`, `generate_asset_tiles.py`, `rgb555()`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 13`** (2 nodes): `dismissActiveHouse()`, `dismissActivePortal()`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `renderFrame()` connect `Community 5` to `Community 0`, `Community 3`, `Community 6`?**
  _High betweenness centrality (0.002) - this node is a cross-community bridge._
- **Why does `clamp()` connect `Community 3` to `Community 0`, `Community 8`, `Community 4`?**
  _High betweenness centrality (0.002) - this node is a cross-community bridge._
- **Why does `init()` connect `Community 6` to `Community 0`, `Community 8`, `Community 5`, `Community 7`?**
  _High betweenness centrality (0.001) - this node is a cross-community bridge._
- **What connects `Screen`, `World`, `House` to the rest of the system?**
  _7 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.07 - nodes in this community are weakly interconnected._