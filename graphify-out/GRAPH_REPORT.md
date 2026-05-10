# Graph Report - fluid-sim  (2026-05-11)

## Corpus Check
- 2 files · ~3,879 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 17 nodes · 31 edges · 2 communities detected
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 2|Community 2]]

## God Nodes (most connected - your core abstractions)
1. `idx()` - 8 edges
2. `project()` - 6 edges
3. `setBoundary()` - 5 edges
4. `linSolve()` - 5 edges
5. `advect()` - 4 edges
6. `step()` - 4 edges
7. `diffuse()` - 3 edges
8. `divergenceAt()` - 3 edges
9. `addDensity()` - 2 edges
10. `addVelocity()` - 2 edges

## Surprising Connections (you probably didn't know these)
- `setBoundary()` --calls--> `idx()`  [EXTRACTED]
  fluid.zig → fluid.zig  _Bridges community 2 → community 0_

## Communities

### Community 0 - "Community 0"
Cohesion: 0.53
Nodes (6): advect(), diffuse(), linSolve(), project(), setBoundary(), step()

### Community 2 - "Community 2"
Cohesion: 0.5
Nodes (4): addDensity(), addVelocity(), divergenceAt(), idx()

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `idx()` connect `Community 2` to `Community 0`, `Community 1`?**
  _High betweenness centrality (0.060) - this node is a cross-community bridge._
- **Why does `project()` connect `Community 0` to `Community 1`, `Community 2`?**
  _High betweenness centrality (0.018) - this node is a cross-community bridge._
- **Why does `linSolve()` connect `Community 0` to `Community 1`, `Community 2`?**
  _High betweenness centrality (0.011) - this node is a cross-community bridge._