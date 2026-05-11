#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from PIL import Image


PALETTE_BASE = 40
TRANSPARENT = 255

TILE_IDS = sorted(
    {
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        60,
        61,
        62,
        63,
        64,
        65,
        66,
        67,
        68,
        72,
        73,
        74,
        75,
        76,
        77,
        78,
        79,
        84,
        85,
        86,
        87,
        88,
        89,
        90,
        91,
        93,
        95,
        96,
        97,
        98,
        99,
        100,
        105,
        107,
        108,
        109,
        110,
        111,
        112,
        120,
        121,
        122,
        123,
        124,
        130,
        131,
    }
)


def rgb555(r: int, g: int, b: int) -> int:
    return ((r >> 3) | ((g >> 3) << 5) | ((b >> 3) << 10))


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    tiles_dir = repo_root / "Asset" / "Tiles"
    if not tiles_dir.exists():
        tiles_dir = repo_root.parent / "Asset" / "Tiles"
    out_file = repo_root / "asset_tiles.zig"

    unique_rgb: dict[tuple[int, int, int], int] = {}
    tile_pixels: list[list[int]] = []

    for tile_id in TILE_IDS:
        tile_path = tiles_dir / f"tile_{tile_id:04d}.png"
        image = Image.open(tile_path).convert("RGBA")
        pixels = list(image.getdata())
        values: list[int] = []
        for r, g, b, a in pixels:
            if a < 128:
                values.append(TRANSPARENT)
                continue

            key = (r, g, b)
            if key not in unique_rgb:
                unique_rgb[key] = len(unique_rgb)
            values.append(PALETTE_BASE + unique_rgb[key])
        tile_pixels.append(values)

    palette_entries = [0] * len(unique_rgb)
    for rgb, idx in unique_rgb.items():
        palette_entries[idx] = rgb555(*rgb)

    lines: list[str] = []
    lines.append("pub const palette_base: u8 = %d;" % PALETTE_BASE)
    lines.append("pub const transparent: u8 = %d;" % TRANSPARENT)
    lines.append("")
    lines.append("pub const tile_ids = [_]u8{")
    for tile_id in TILE_IDS:
        lines.append(f"    {tile_id},")
    lines.append("};")
    lines.append("")
    lines.append("pub const palette = [_]u16{")
    for value in palette_entries:
        lines.append(f"    0x{value:04x},")
    lines.append("};")
    lines.append("")
    lines.append("pub const tiles = [_][16 * 16]u8{")
    for tile in tile_pixels:
        lines.append("    .{")
        for row in range(16):
            start = row * 16
            row_values = ", ".join(str(v) for v in tile[start : start + 16])
            lines.append(f"        {row_values},")
        lines.append("    },")
    lines.append("};")
    lines.append("")

    out_file.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_file}")
    print(f"tile_count={len(TILE_IDS)} palette_entries={len(palette_entries)}")


if __name__ == "__main__":
    main()
