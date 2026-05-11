const std = @import("std");
const asset = @import("asset_tiles.zig");

const Screen = struct {
    const width: usize = 240;
    const height: usize = 160;
};

const World = struct {
    const tile_size: i32 = 16;
    const map_w: usize = 64;
    const map_h: usize = 44;
    const width: f32 = @as(f32, @floatFromInt(map_w * @as(usize, tile_size)));
    const height: f32 = @as(f32, @floatFromInt(map_h * @as(usize, tile_size)));
};

const N: usize = 48;
const SIZE: usize = (N + 2) * (N + 2);
const ITER: usize = 4;
const HOUSE_COUNT: usize = 5;
const RIPPLE_COUNT: usize = 72;
const PARTICLE_COUNT: usize = 72;
const FRAMEBUFFER_SIZE: usize = Screen.width * Screen.height;
const PALETTE_SIZE: usize = 256;
const MAP_TILE_COUNT: usize = World.map_w * World.map_h;

const BUTTON_ACCEPT: u32 = 1 << 0;
const BUTTON_CANCEL: u32 = 1 << 1;
const BUTTON_REDUCED_MOTION: u32 = 1 << 2;

const PLAYER_SPEED: f32 = 74.0;
const TRIGGER_RADIUS: f32 = 78.0;
const HYSTERESIS: f32 = 24.0;
const DIFFUSION: f32 = 0.001;
const VISCOSITY: f32 = 0.0001;
const DT: f32 = 0.2;
const FADE: f32 = 0.992;
const MIN_DENSITY: f32 = 0.01;
const MAX_DENSITY: f32 = 1000.0;

const House = struct {
    x: f32,
    y: f32,
    radius: f32,
};

const Player = struct {
    x: f32,
    y: f32,
    vx: f32,
    vy: f32,
    dir: f32,
    moving: f32,
    step_clock: f32,
    ripple_clock: f32,
};

const Ripple = struct {
    active: bool,
    x: f32,
    y: f32,
    radius: f32,
    life: f32,
    alpha: f32,
};

const Particle = struct {
    active: bool,
    x: f32,
    y: f32,
    vx: f32,
    vy: f32,
    life: f32,
};

const houses = [_]House{
    .{ .x = 224.0, .y = 200.0, .radius = TRIGGER_RADIUS },
    .{ .x = 416.0, .y = 200.0, .radius = TRIGGER_RADIUS },
    .{ .x = 608.0, .y = 200.0, .radius = TRIGGER_RADIUS },
    .{ .x = 320.0, .y = 408.0, .radius = TRIGGER_RADIUS },
    .{ .x = 544.0, .y = 408.0, .radius = TRIGGER_RADIUS },
};

const routes = [_]usize{ 0, 1, 2, 4, 3, 0 };

const portal_brand_bg = [_]u8{ 200, 202, 204, 206, 209 };
const portal_brand_fg = [_]u8{ 201, 201, 205, 208, 201 };

const HouseSprite = struct {
    cols: usize,
    rows: usize,
    tiles: []const u8,
    base_w: i32,
};

const house_sprites = [_]HouseSprite{
    .{
        .cols = 4,
        .rows = 4,
        .tiles = &[_]u8{ 52, 53, 54, 55, 64, 65, 66, 67, 72, 73, 74, 75, 84, 85, 86, 87 },
        .base_w = 44,
    },
    .{
        .cols = 4,
        .rows = 4,
        .tiles = &[_]u8{ 48, 49, 50, 51, 60, 61, 62, 63, 76, 77, 78, 79, 88, 89, 90, 91 },
        .base_w = 44,
    },
    .{
        .cols = 3,
        .rows = 4,
        .tiles = &[_]u8{ 52, 53, 54, 64, 67, 68, 72, 74, 75, 84, 85, 87 },
        .base_w = 38,
    },
    .{
        .cols = 3,
        .rows = 4,
        .tiles = &[_]u8{ 48, 49, 50, 60, 63, 62, 76, 77, 79, 88, 89, 91 },
        .base_w = 38,
    },
    .{
        .cols = 5,
        .rows = 3,
        .tiles = &[_]u8{ 96, 97, 98, 99, 100, 108, 109, 110, 111, 112, 120, 121, 122, 123, 124 },
        .base_w = 56,
    },
};

var player: Player = undefined;
var camera: [2]f32 = .{ 0.0, 0.0 };
var active_house: i32 = -1;
var dismissed_house: i32 = -1;

var density: [SIZE]f32 = undefined;
var density_prev: [SIZE]f32 = undefined;
var vx: [SIZE]f32 = undefined;
var vy: [SIZE]f32 = undefined;
var vx_prev: [SIZE]f32 = undefined;
var vy_prev: [SIZE]f32 = undefined;
var total_density: f32 = 0;

var ripples: [RIPPLE_COUNT]Ripple = undefined;
var particles: [PARTICLE_COUNT]Particle = undefined;

var house_view: [HOUSE_COUNT * 3]f32 = undefined;
var player_view: [6]f32 = undefined;
var ripple_view: [RIPPLE_COUNT * 5]f32 = undefined;
var particle_view: [PARTICLE_COUNT * 5]f32 = undefined;

var palette: [PALETTE_SIZE]u16 = [_]u16{0} ** PALETTE_SIZE;
var framebuffer: [FRAMEBUFFER_SIZE]u16 = [_]u16{0} ** FRAMEBUFFER_SIZE;
var world_tiles: [MAP_TILE_COUNT]u8 = [_]u8{0} ** MAP_TILE_COUNT;

var input_axis_x: f32 = 0.0;
var input_axis_y: f32 = 0.0;
var input_buttons: u32 = 0;
var render_clock: f32 = 0.0;

fn clamp(value: f32, min: f32, max: f32) f32 {
    return @max(min, @min(max, value));
}

fn clampI(value: i32, min: i32, max: i32) i32 {
    return @max(min, @min(max, value));
}

fn idx(x: usize, y: usize) usize {
    return x + (N + 2) * y;
}

fn distance(ax: f32, ay: f32, bx: f32, by: f32) f32 {
    const dx = ax - bx;
    const dy = ay - by;
    return @sqrt(dx * dx + dy * dy);
}

fn rgb555(r: u8, g: u8, b: u8) u16 {
    const rr: u16 = @as(u16, r) >> 3;
    const gg: u16 = @as(u16, g) >> 3;
    const bb: u16 = @as(u16, b) >> 3;
    return rr | (gg << 5) | (bb << 10);
}

fn initPalette() void {
    palette[0] = rgb555(20, 40, 52);
    palette[1] = rgb555(23, 78, 110);
    palette[2] = rgb555(30, 98, 138);
    palette[3] = rgb555(43, 118, 162);
    palette[4] = rgb555(58, 142, 180);
    palette[5] = rgb555(99, 176, 201);
    palette[6] = rgb555(128, 206, 217);
    palette[7] = rgb555(160, 231, 225);
    palette[8] = rgb555(16, 32, 40);
    palette[9] = rgb555(48, 70, 86);
    palette[10] = rgb555(92, 55, 37);
    palette[11] = rgb555(132, 79, 51);
    palette[12] = rgb555(169, 113, 69);
    palette[13] = rgb555(214, 156, 94);
    palette[14] = rgb555(56, 62, 73);
    palette[15] = rgb555(128, 136, 152);
    palette[16] = rgb555(186, 194, 205);
    palette[17] = rgb555(204, 170, 112);
    palette[18] = rgb555(238, 202, 131);
    palette[19] = rgb555(238, 207, 165);
    palette[20] = rgb555(132, 60, 78);
    palette[21] = rgb555(255, 244, 208);
    palette[22] = rgb555(226, 115, 145);
    palette[23] = rgb555(74, 92, 109);
    palette[24] = rgb555(177, 86, 65);
    palette[25] = rgb555(211, 115, 82);
    palette[26] = rgb555(62, 54, 61);

    // Brand hues for portal pinpoints.
    palette[200] = rgb555(10, 102, 194); // LinkedIn
    palette[201] = rgb555(247, 250, 252); // bright white
    palette[202] = rgb555(24, 23, 23); // GitHub
    palette[203] = rgb555(255, 0, 80); // TikTok pink
    palette[204] = rgb555(228, 64, 95); // Instagram pink
    palette[205] = rgb555(252, 175, 69); // Instagram warm accent
    palette[206] = rgb555(0, 242, 234); // TikTok cyan
    palette[207] = rgb555(67, 97, 238); // Portfolio blue
    palette[208] = rgb555(246, 190, 96); // Portfolio amber
    palette[209] = rgb555(58, 12, 163); // Instagram deep violet

    for (asset.palette, 0..) |color, i| {
        const slot = @as(usize, asset.palette_base) + i;
        if (slot < PALETTE_SIZE) {
            palette[slot] = color;
        }
    }
}

fn findAssetTile(sheet_id: u8) ?usize {
    for (asset.tile_ids, 0..) |tile_id, i| {
        if (tile_id == sheet_id) return i;
    }
    return null;
}

fn drawAssetTile16(sheet_id: u8, x: i32, y: i32) void {
    const tile_idx = findAssetTile(sheet_id) orelse return;
    const tile = asset.tiles[tile_idx];

    var sy: usize = 0;
    while (sy < 16) : (sy += 1) {
        var sx: usize = 0;
        while (sx < 16) : (sx += 1) {
            const color = tile[sy * 16 + sx];
            if (color == asset.transparent) continue;
            setPixel(x + @as(i32, @intCast(sx)), y + @as(i32, @intCast(sy)), color);
        }
    }
}

fn mapOffset(tx: i32, ty: i32) ?usize {
    const map_w_i: i32 = @intCast(World.map_w);
    const map_h_i: i32 = @intCast(World.map_h);
    if (tx < 0 or ty < 0) return null;
    if (tx >= map_w_i or ty >= map_h_i) return null;
    return @as(usize, @intCast(tx)) + @as(usize, @intCast(ty)) * World.map_w;
}

fn setMapTile(tx: i32, ty: i32, tile_id: u8) void {
    const offset = mapOffset(tx, ty) orelse return;
    world_tiles[offset] = tile_id;
}

fn fillMapRect(tx: i32, ty: i32, tw: i32, th: i32, tile_id: u8) void {
    var y = ty;
    while (y < ty + th) : (y += 1) {
        var x = tx;
        while (x < tx + tw) : (x += 1) {
            setMapTile(x, y, tile_id);
        }
    }
}

fn drawMapTile(tx: i32, ty: i32) void {
    const offset = mapOffset(tx, ty) orelse return;
    const tile_id = world_tiles[offset];
    const sx = tx * World.tile_size - @as(i32, @intFromFloat(camera[0]));
    const sy = ty * World.tile_size - @as(i32, @intFromFloat(camera[1]));
    drawAssetTile16(tile_id, sx, sy);
}

fn initTownMap() void {
    for (0..World.map_h) |ty| {
        for (0..World.map_w) |tx| {
            const jitter = (tx * 7 + ty * 11) % 9;
            const tile_id: u8 = if (jitter == 0) 2 else if (jitter < 3) 1 else 0;
            setMapTile(@as(i32, @intCast(tx)), @as(i32, @intCast(ty)), tile_id);
        }
    }

    // Main streets.
    fillMapRect(7, 11, 50, 3, 13);
    fillMapRect(11, 23, 44, 3, 13);
    fillMapRect(13, 7, 3, 27, 13);
    fillMapRect(25, 7, 3, 27, 13);
    fillMapRect(37, 7, 3, 27, 13);

    // Side blocks for neighborhoods.
    fillMapRect(6, 16, 10, 5, 12);
    fillMapRect(22, 16, 10, 5, 12);
    fillMapRect(38, 16, 10, 5, 12);
    fillMapRect(12, 28, 10, 5, 12);
    fillMapRect(30, 28, 10, 5, 12);

    // Fences and props around plots.
    var x: i32 = 5;
    while (x < 58) : (x += 4) {
        setMapTile(x, 10, 69);
        setMapTile(x, 26, 69);
    }
    setMapTile(5, 10, 68);
    setMapTile(57, 10, 70);
    setMapTile(5, 26, 68);
    setMapTile(57, 26, 70);

    // Tree rows for Medabots-like suburban blocks.
    for (0..7) |i| {
        setMapTile(3 + @as(i32, @intCast(i * 2)), 5, if (i % 2 == 0) 4 else 3);
        setMapTile(47 + @as(i32, @intCast(i * 2)), 5, if (i % 2 == 0) 10 else 11);
    }
    for (0..10) |i| {
        setMapTile(2 + @as(i32, @intCast(i * 2)), 36, if (i % 3 == 0) 5 else 9);
    }

    // Utility objects and signs around the streets.
    setMapTile(20, 14, 82);
    setMapTile(35, 14, 82);
    setMapTile(23, 24, 94);
    setMapTile(41, 24, 94);
    setMapTile(17, 30, 83);
    setMapTile(43, 30, 83);
}

fn drawTownMap() void {
    const tile_size = World.tile_size;
    const map_w_i: i32 = @intCast(World.map_w);
    const map_h_i: i32 = @intCast(World.map_h);
    const visible_tiles_x: i32 = @as(i32, @intCast(Screen.width / 16));
    const visible_tiles_y: i32 = @as(i32, @intCast(Screen.height / 16));

    const start_tx = clampI(@divFloor(@as(i32, @intFromFloat(camera[0])), tile_size) - 1, 0, map_w_i - 1);
    const start_ty = clampI(@divFloor(@as(i32, @intFromFloat(camera[1])), tile_size) - 1, 0, map_h_i - 1);
    const end_tx = clampI(start_tx + visible_tiles_x + 3, 0, map_w_i - 1);
    const end_ty = clampI(start_ty + visible_tiles_y + 3, 0, map_h_i - 1);

    var ty = start_ty;
    while (ty <= end_ty) : (ty += 1) {
        var tx = start_tx;
        while (tx <= end_tx) : (tx += 1) {
            drawMapTile(tx, ty);
        }
    }
}

fn drawRouteDots() void {
    var route_i: usize = 1;
    while (route_i < routes.len) : (route_i += 1) {
        const from = houses[routes[route_i - 1]];
        const to = houses[routes[route_i]];
        const dx = to.x - from.x;
        const dy = to.y - from.y;
        const distance_len = @sqrt(dx * dx + dy * dy);
        const steps = @max(@as(usize, 1), @as(usize, @intFromFloat(distance_len / 32.0)));

        var step: usize = 1;
        while (step <= steps) : (step += 1) {
            const t = @as(f32, @floatFromInt(step)) / @as(f32, @floatFromInt(steps + 1));
            const world_x = from.x + dx * t;
            const world_y = from.y + dy * t;
            const sx = @as(i32, @intFromFloat(world_x - camera[0]));
            const sy = @as(i32, @intFromFloat(world_y - camera[1]));
            fillRect(sx - 1, sy - 1, 3, 3, 21);
            fillRect(sx, sy, 1, 1, 8);
        }
    }
}

fn drawPortalBadge(index: usize, anchor_x: i32, anchor_y: i32) void {
    const bg = portal_brand_bg[index % portal_brand_bg.len];
    const fg = portal_brand_fg[index % portal_brand_fg.len];

    fillRect(anchor_x - 7, anchor_y - 9, 14, 10, 8);
    fillRect(anchor_x - 6, anchor_y - 8, 12, 8, bg);

    switch (index) {
        // LinkedIn style `in`
        0 => {
            fillRect(anchor_x - 4, anchor_y - 6, 2, 6, fg);
            fillRect(anchor_x - 4, anchor_y - 8, 2, 1, fg);
            fillRect(anchor_x - 1, anchor_y - 4, 2, 4, fg);
            fillRect(anchor_x + 1, anchor_y - 6, 2, 6, fg);
            fillRect(anchor_x + 1, anchor_y - 5, 2, 1, fg);
        },
        // GitHub-like cat face silhouette
        1 => {
            fillRect(anchor_x - 3, anchor_y - 6, 8, 6, fg);
            fillRect(anchor_x - 4, anchor_y - 5, 1, 4, fg);
            fillRect(anchor_x + 4, anchor_y - 5, 1, 4, fg);
            fillRect(anchor_x - 2, anchor_y - 7, 2, 1, fg);
            fillRect(anchor_x + 1, anchor_y - 7, 2, 1, fg);
            fillRect(anchor_x - 1, anchor_y - 3, 1, 1, bg);
            fillRect(anchor_x + 1, anchor_y - 3, 1, 1, bg);
        },
        // TikTok note silhouette
        2 => {
            fillRect(anchor_x + 1, anchor_y - 7, 2, 6, fg);
            fillRect(anchor_x - 2, anchor_y - 4, 5, 2, fg);
            fillRect(anchor_x - 3, anchor_y - 2, 4, 2, fg);
        },
        // Instagram camera glyph
        3 => {
            fillRect(anchor_x - 4, anchor_y - 7, 8, 1, fg);
            fillRect(anchor_x - 4, anchor_y - 1, 8, 1, fg);
            fillRect(anchor_x - 4, anchor_y - 6, 1, 5, fg);
            fillRect(anchor_x + 3, anchor_y - 6, 1, 5, fg);
            fillRect(anchor_x - 1, anchor_y - 5, 2, 2, 205);
            fillRect(anchor_x + 1, anchor_y - 6, 1, 1, fg);
        },
        // Portfolio briefcase icon
        else => {
            fillRect(anchor_x - 4, anchor_y - 5, 8, 5, fg);
            fillRect(anchor_x - 2, anchor_y - 7, 4, 2, fg);
            fillRect(anchor_x - 1, anchor_y - 4, 2, 1, 208);
        },
    }
}

fn drawHouse(index: usize) void {
    const house = houses[index];
    const sx = @as(i32, @intFromFloat(house.x - camera[0]));
    const sy = @as(i32, @intFromFloat(house.y - camera[1]));
    const sprite = house_sprites[index % house_sprites.len];

    const house_w: i32 = @as(i32, @intCast(sprite.cols)) * 16;
    const house_h: i32 = @as(i32, @intCast(sprite.rows)) * 16;
    const house_x = sx - @divTrunc(house_w, 2);
    const house_y = sy - house_h + 10;

    fillRect(house_x + 4, sy + 2, house_w - 8, 6, 8);

    var row: usize = 0;
    while (row < sprite.rows) : (row += 1) {
        var col: usize = 0;
        while (col < sprite.cols) : (col += 1) {
            const tile_id = sprite.tiles[row * sprite.cols + col];
            drawAssetTile16(tile_id, house_x + @as(i32, @intCast(col * 16)), house_y + @as(i32, @intCast(row * 16)));
        }
    }

    drawPortalBadge(index, sx, house_y - 4);

    if (active_house == @as(i32, @intCast(index))) {
        drawFrameRect(house_x - 2, house_y - 2, house_w + 4, house_h + 6, 21);
    }
}

fn syncHouseView() void {
    for (houses, 0..) |house, i| {
        const base = i * 3;
        house_view[base] = house.x;
        house_view[base + 1] = house.y;
        house_view[base + 2] = house.radius;
    }
}

fn syncPlayerView() void {
    player_view[0] = player.x;
    player_view[1] = player.y;
    player_view[2] = player.vx;
    player_view[3] = player.vy;
    player_view[4] = player.dir;
    player_view[5] = player.moving;
}

fn syncRippleView() void {
    for (ripples, 0..) |ripple, i| {
        const base = i * 5;
        ripple_view[base] = ripple.x;
        ripple_view[base + 1] = ripple.y;
        ripple_view[base + 2] = ripple.radius;
        ripple_view[base + 3] = ripple.alpha;
        ripple_view[base + 4] = if (ripple.active) 1.0 else 0.0;
    }
}

fn syncParticleView() void {
    for (particles, 0..) |particle, i| {
        const base = i * 5;
        particle_view[base] = particle.x;
        particle_view[base + 1] = particle.y;
        particle_view[base + 2] = particle.vx;
        particle_view[base + 3] = particle.vy;
        particle_view[base + 4] = if (particle.active) particle.life else 0.0;
    }
}

fn syncViews() void {
    syncPlayerView();
    syncRippleView();
    syncParticleView();
}

fn resetFluid() void {
    @memset(&density, 0);
    @memset(&density_prev, 0);
    @memset(&vx, 0);
    @memset(&vy, 0);
    @memset(&vx_prev, 0);
    @memset(&vy_prev, 0);
    total_density = 0;
}

fn resetEffects() void {
    for (&ripples) |*ripple| {
        ripple.* = .{ .active = false, .x = 0, .y = 0, .radius = 0, .life = 0, .alpha = 0 };
    }
    for (&particles) |*particle| {
        particle.* = .{ .active = false, .x = 0, .y = 0, .vx = 0, .vy = 0, .life = 0 };
    }
}

fn reducedMotionEnabled() bool {
    return (input_buttons & BUTTON_REDUCED_MOTION) != 0;
}

fn fbOffset(x: i32, y: i32) ?usize {
    if (x < 0 or y < 0) return null;
    if (x >= Screen.width or y >= Screen.height) return null;
    return @as(usize, @intCast(x)) + @as(usize, @intCast(y)) * Screen.width;
}

fn setPixel(x: i32, y: i32, palette_index: u8) void {
    const offset = fbOffset(x, y) orelse return;
    framebuffer[offset] = palette[palette_index];
}

fn fillRect(x: i32, y: i32, w: i32, h: i32, palette_index: u8) void {
    if (w <= 0 or h <= 0) return;
    const min_x = clampI(x, 0, Screen.width);
    const min_y = clampI(y, 0, Screen.height);
    const max_x = clampI(x + w, 0, Screen.width);
    const max_y = clampI(y + h, 0, Screen.height);
    if (min_x >= max_x or min_y >= max_y) return;

    var py = min_y;
    while (py < max_y) : (py += 1) {
        const row_base = @as(usize, @intCast(py)) * Screen.width;
        var px = min_x;
        while (px < max_x) : (px += 1) {
            const offset = row_base + @as(usize, @intCast(px));
            framebuffer[offset] = palette[palette_index];
        }
    }
}

fn drawFrameRect(x: i32, y: i32, w: i32, h: i32, palette_index: u8) void {
    fillRect(x, y, w, 1, palette_index);
    fillRect(x, y + h - 1, w, 1, palette_index);
    fillRect(x, y, 1, h, palette_index);
    fillRect(x + w - 1, y, 1, h, palette_index);
}

fn drawPlayer() void {
    const sx = @as(i32, @intFromFloat(player.x - camera[0]));
    const sy = @as(i32, @intFromFloat(player.y - camera[1]));
    const moving = player.moving > 0.5;

    const frame = if (moving) @as(i32, @intFromFloat(@floor(player.step_clock * 2.0))) & 1 else 0;
    const leg_offset: i32 = if (frame == 0) 1 else -1;

    fillRect(sx - 7, sy + 2, 14, 3, 8);
    fillRect(sx - 4, sy - 16, 8, 6, 19);
    fillRect(sx - 6, sy - 10, 12, 8, 24);
    fillRect(sx - 5, sy - 2, 4, 6, 18);
    fillRect(sx + 1, sy - 2 + leg_offset, 4, 6, 18);
    fillRect(sx - 5, sy + 4, 4, 1, 8);
    fillRect(sx + 1, sy + 4 + leg_offset, 4, 1, 8);

    if (player.dir == 2.0) {
        fillRect(sx - 3, sy - 14, 2, 2, 8);
    } else if (player.dir == 3.0) {
        fillRect(sx + 1, sy - 14, 2, 2, 8);
    } else {
        fillRect(sx - 2, sy - 14, 2, 2, 8);
        fillRect(sx + 1, sy - 14, 2, 2, 8);
    }
}

fn drawMiniMap() void {
    const map_w: i32 = 66;
    const map_h: i32 = 44;
    const map_x: i32 = @as(i32, @intCast(Screen.width)) - map_w - 6;
    const map_y: i32 = @as(i32, @intCast(Screen.height)) - map_h - 6;

    fillRect(map_x, map_y, map_w, map_h, 8);
    fillRect(map_x + 1, map_y + 1, map_w - 2, map_h - 2, 0);
    drawFrameRect(map_x, map_y, map_w, map_h, 21);

    var i: usize = 0;
    while (i < HOUSE_COUNT) : (i += 1) {
        const hx = map_x + 2 + @as(i32, @intFromFloat((houses[i].x / World.width) * @as(f32, @floatFromInt(map_w - 4))));
        const hy = map_y + 2 + @as(i32, @intFromFloat((houses[i].y / World.height) * @as(f32, @floatFromInt(map_h - 4))));
        fillRect(hx - 1, hy - 1, 4, 4, portal_brand_bg[i % portal_brand_bg.len]);
        fillRect(hx, hy, 2, 2, portal_brand_fg[i % portal_brand_fg.len]);
    }

    const px = map_x + 2 + @as(i32, @intFromFloat((player.x / World.width) * @as(f32, @floatFromInt(map_w - 4))));
    const py = map_y + 2 + @as(i32, @intFromFloat((player.y / World.height) * @as(f32, @floatFromInt(map_h - 4))));
    fillRect(px - 1, py - 1, 4, 4, 8);
    fillRect(px, py, 2, 2, 21);

    const view_x = map_x + 2 + @as(i32, @intFromFloat((camera[0] / World.width) * @as(f32, @floatFromInt(map_w - 4))));
    const view_y = map_y + 2 + @as(i32, @intFromFloat((camera[1] / World.height) * @as(f32, @floatFromInt(map_h - 4))));
    const view_w = @max(4, @as(i32, @intFromFloat((@as(f32, @floatFromInt(Screen.width)) / World.width) * @as(f32, @floatFromInt(map_w - 4)))));
    const view_h = @max(4, @as(i32, @intFromFloat((@as(f32, @floatFromInt(Screen.height)) / World.height) * @as(f32, @floatFromInt(map_h - 4)))));
    drawFrameRect(view_x, view_y, view_w, view_h, 6);
}

fn renderFrame() void {
    fillRect(0, 0, Screen.width, Screen.height, 0);

    drawTownMap();

    drawRouteDots();

    var draw_list: [HOUSE_COUNT + 1]struct { is_player: bool, y: f32, index: usize } = undefined;
    draw_list[0] = .{ .is_player = true, .y = player.y, .index = 0 };
    var i: usize = 0;
    while (i < HOUSE_COUNT) : (i += 1) {
        draw_list[i + 1] = .{ .is_player = false, .y = houses[i].y, .index = i };
    }

    std.mem.sort(@TypeOf(draw_list[0]), &draw_list, {}, struct {
        fn lessThan(_: void, a: @TypeOf(draw_list[0]), b: @TypeOf(draw_list[0])) bool {
            return a.y < b.y;
        }
    }.lessThan);

    for (draw_list) |item| {
        if (item.is_player) {
            drawPlayer();
        } else {
            drawHouse(item.index);
        }
    }

    drawMiniMap();

    fillRect(4, 4, 44, 12, 8);
    drawFrameRect(4, 4, 44, 12, 21);
    fillRect(5, 5, 42, 10, 23);

    fillRect(@as(i32, @intCast(Screen.width)) - 46, 4, 42, 12, 8);
    drawFrameRect(@as(i32, @intCast(Screen.width)) - 46, 4, 42, 12, 21);
}

export fn init() void {
    player = .{
        .x = 416.0,
        .y = 256.0,
        .vx = 0.0,
        .vy = 0.0,
        .dir = 0.0,
        .moving = 0.0,
        .step_clock = 0.0,
        .ripple_clock = 0.0,
    };
    camera = .{ 296.0, 176.0 };
    active_house = -1;
    dismissed_house = -1;
    input_axis_x = 0.0;
    input_axis_y = 0.0;
    input_buttons = 0;
    render_clock = 0.0;
    initTownMap();
    resetFluid();
    resetEffects();
    initPalette();
    syncHouseView();
    syncViews();
    renderFrame();
}

export fn getHouseCount() usize {
    return HOUSE_COUNT;
}

export fn getRippleCount() usize {
    return RIPPLE_COUNT;
}

export fn getParticleCount() usize {
    return PARTICLE_COUNT;
}

export fn getHousePtr() [*]u8 {
    return @ptrCast(&house_view);
}

export fn getPlayerPtr() [*]u8 {
    return @ptrCast(&player_view);
}

export fn getCameraPtr() [*]u8 {
    return @ptrCast(&camera);
}

export fn getRipplePtr() [*]u8 {
    return @ptrCast(&ripple_view);
}

export fn getParticlePtr() [*]u8 {
    return @ptrCast(&particle_view);
}

export fn getDensityPtr() [*]u8 {
    return @ptrCast(&density);
}

export fn getGridSize() usize {
    return N;
}

export fn getWorldWidth() f32 {
    return World.width;
}

export fn getWorldHeight() f32 {
    return World.height;
}

export fn getTotalDensity() f32 {
    return total_density;
}

export fn getFramebufferPtr() [*]u16 {
    return @ptrCast(&framebuffer);
}

export fn getPalettePtr() [*]u16 {
    return @ptrCast(&palette);
}

export fn getFramebufferLen() usize {
    return FRAMEBUFFER_SIZE;
}

export fn getActivePortal() i32 {
    return active_house;
}

export fn getActiveHouse() i32 {
    return active_house;
}

export fn dismissActivePortal() void {
    dismissed_house = active_house;
    active_house = -1;
}

export fn dismissActiveHouse() void {
    dismissActivePortal();
}

export fn setInput(button_mask: u32, axis_x: f32, axis_y: f32) void {
    input_buttons = button_mask;
    input_axis_x = clamp(axis_x, -1.0, 1.0);
    input_axis_y = clamp(axis_y, -1.0, 1.0);
}

export fn setPlayerPosition(x: f32, y: f32) void {
    player.x = clamp(x, 26.0, World.width - 26.0);
    player.y = clamp(y, 34.0, World.height - 26.0);
    syncPlayerView();
}

export fn addDensity(x: usize, y: usize, amount: f32) void {
    if (x >= 1 and x <= N and y >= 1 and y <= N) {
        density[idx(x, y)] += amount;
    }
}

export fn addVelocity(x: usize, y: usize, amount_x: f32, amount_y: f32) void {
    if (x >= 1 and x <= N and y >= 1 and y <= N) {
        vx[idx(x, y)] += amount_x;
        vy[idx(x, y)] += amount_y;
    }
}

fn spawnRipple(x: f32, y: f32, force: f32, reduced_motion: bool) void {
    if (reduced_motion) return;
    for (&ripples) |*ripple| {
        if (!ripple.active) {
            ripple.* = .{
                .active = true,
                .x = x,
                .y = y + 8.0,
                .radius = 3.0,
                .life = 1.0,
                .alpha = clamp(0.32 + force * 0.18, 0.28, 0.56),
            };
            return;
        }
    }
}

fn spawnSplash(x: f32, y: f32, reduced_motion: bool) void {
    if (reduced_motion) return;
    var made: usize = 0;
    for (&particles) |*particle| {
        if (particle.active) continue;
        const angle = (@as(f32, @floatFromInt(made)) / 6.0) * std.math.tau;
        particle.* = .{
            .active = true,
            .x = x,
            .y = y,
            .vx = @cos(angle) * 22.0,
            .vy = @sin(angle) * 12.0 - 12.0,
            .life = 0.34,
        };
        made += 1;
        if (made == 6) return;
    }
}

fn disturbFluid(world_x: f32, world_y: f32, amount_x: f32, amount_y: f32) void {
    const gx_float = clamp((world_x / World.width) * @as(f32, @floatFromInt(N)), 1.0, @as(f32, @floatFromInt(N)));
    const gy_float = clamp((world_y / World.height) * @as(f32, @floatFromInt(N)), 1.0, @as(f32, @floatFromInt(N)));
    const gx: usize = @intFromFloat(gx_float);
    const gy: usize = @intFromFloat(gy_float);
    addDensity(gx, gy, 76.0);
    addVelocity(gx, gy, amount_x * 0.08, amount_y * 0.08);
}

fn updatePlayer(dt: f32, input_x_raw: f32, input_y_raw: f32, reduced_motion: bool) void {
    _ = reduced_motion;
    var input_x = clamp(input_x_raw, -1.0, 1.0);
    var input_y = clamp(input_y_raw, -1.0, 1.0);
    const input_len = distance(input_x, input_y, 0.0, 0.0);
    if (input_len > 1.0) {
        input_x /= input_len;
        input_y /= input_len;
    }

    player.vx = input_x * PLAYER_SPEED;
    player.vy = input_y * PLAYER_SPEED;
    const moving = input_len > 0.05;
    player.moving = if (moving) 1.0 else 0.0;

    if (@abs(input_x) > @abs(input_y)) {
        player.dir = if (input_x < 0) 2.0 else 3.0;
    } else if (@abs(input_y) > 0.05) {
        player.dir = if (input_y < 0) 1.0 else 0.0;
    }

    player.x = clamp(player.x + player.vx * dt, 26.0, World.width - 26.0);
    player.y = clamp(player.y + player.vy * dt, 34.0, World.height - 26.0);
    const step_rate: f32 = if (moving) 9.0 else 2.0;
    player.step_clock += dt * step_rate;
}

fn updateCamera() void {
    const target_x = clamp(player.x - @as(f32, @floatFromInt(Screen.width)) / 2.0, 0.0, World.width - @as(f32, @floatFromInt(Screen.width)));
    const target_y = clamp(player.y - @as(f32, @floatFromInt(Screen.height)) / 2.0, 0.0, World.height - @as(f32, @floatFromInt(Screen.height)));
    camera[0] += (target_x - camera[0]) * 0.16;
    camera[1] += (target_y - camera[1]) * 0.16;
}

fn updateEffects(dt: f32) void {
    for (&ripples) |*ripple| {
        if (!ripple.active) continue;
        ripple.radius += 38.0 * dt;
        ripple.life -= 1.7 * dt;
        ripple.alpha *= 1.0 - 2.6 * dt;
        if (ripple.life <= 0.0 or ripple.alpha <= 0.02) ripple.active = false;
    }

    for (&particles) |*particle| {
        if (!particle.active) continue;
        particle.life -= dt;
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
        particle.vy += 44.0 * dt;
        if (particle.life <= 0.0) particle.active = false;
    }
}

fn updateProximity() void {
    if (active_house >= 0) {
        const house = houses[@intCast(active_house)];
        if (distance(player.x, player.y, house.x, house.y) > TRIGGER_RADIUS + HYSTERESIS) {
            active_house = -1;
        }
        return;
    }

    var nearest: i32 = -1;
    var nearest_distance: f32 = 99999.0;
    for (houses, 0..) |house, i| {
        const d = distance(player.x, player.y, house.x, house.y);
        if (d < nearest_distance) {
            nearest_distance = d;
            nearest = @intCast(i);
        }
    }

    if (dismissed_house >= 0) {
        const dismissed = houses[@intCast(dismissed_house)];
        if (distance(player.x, player.y, dismissed.x, dismissed.y) > TRIGGER_RADIUS + HYSTERESIS) {
            dismissed_house = -1;
        }
    }

    if (nearest >= 0 and nearest_distance < TRIGGER_RADIUS and nearest != dismissed_house) {
        active_house = nearest;
    }
}

fn setBoundary(b: i32, x: []f32) void {
    var i: usize = 1;
    while (i <= N) : (i += 1) {
        x[idx(0, i)] = if (b == 1) -x[idx(1, i)] else x[idx(1, i)];
        x[idx(N + 1, i)] = if (b == 1) -x[idx(N, i)] else x[idx(N, i)];
        x[idx(i, 0)] = if (b == 2) -x[idx(i, 1)] else x[idx(i, 1)];
        x[idx(i, N + 1)] = if (b == 2) -x[idx(i, N)] else x[idx(i, N)];
    }

    x[idx(0, 0)] = 0.5 * (x[idx(1, 0)] + x[idx(0, 1)]);
    x[idx(0, N + 1)] = 0.5 * (x[idx(1, N + 1)] + x[idx(0, N)]);
    x[idx(N + 1, 0)] = 0.5 * (x[idx(N, 0)] + x[idx(N + 1, 1)]);
    x[idx(N + 1, N + 1)] = 0.5 * (x[idx(N, N + 1)] + x[idx(N + 1, N)]);
}

fn linSolve(b: i32, x: []f32, x0: []f32, a: f32, c: f32) void {
    var k: usize = 0;
    while (k < ITER) : (k += 1) {
        var j: usize = 1;
        while (j <= N) : (j += 1) {
            var i: usize = 1;
            while (i <= N) : (i += 1) {
                x[idx(i, j)] = (x0[idx(i, j)] + a * (x[idx(i - 1, j)] + x[idx(i + 1, j)] + x[idx(i, j - 1)] + x[idx(i, j + 1)])) / c;
            }
        }
        setBoundary(b, x);
    }
}

fn diffuse(b: i32, x: []f32, x0: []f32, diff: f32, dt: f32) void {
    const a = dt * diff * @as(f32, @floatFromInt(N * N));
    linSolve(b, x, x0, a, 1.0 + 4.0 * a);
}

fn divergenceAt(u: []const f32, v: []const f32, i: usize, j: usize) f32 {
    const inv_grid = 1.0 / @as(f32, @floatFromInt(N));
    return -0.5 * inv_grid * (u[idx(i + 1, j)] - u[idx(i - 1, j)] + v[idx(i, j + 1)] - v[idx(i, j - 1)]);
}

fn advect(b: i32, d: []f32, d0: []f32, u: []f32, v: []f32, dt: f32) void {
    const dt0 = dt * @as(f32, @floatFromInt(N));

    var j: usize = 1;
    while (j <= N) : (j += 1) {
        var i: usize = 1;
        while (i <= N) : (i += 1) {
            var x = @as(f32, @floatFromInt(i)) - dt0 * u[idx(i, j)];
            var y = @as(f32, @floatFromInt(j)) - dt0 * v[idx(i, j)];

            x = clamp(x, 0.5, @as(f32, @floatFromInt(N)) + 0.5);
            y = clamp(y, 0.5, @as(f32, @floatFromInt(N)) + 0.5);

            const ix0: usize = @intFromFloat(x);
            const ix1 = ix0 + 1;
            const iy0: usize = @intFromFloat(y);
            const iy1 = iy0 + 1;

            const s1 = x - @as(f32, @floatFromInt(ix0));
            const s0 = 1.0 - s1;
            const t1 = y - @as(f32, @floatFromInt(iy0));
            const t0 = 1.0 - t1;

            d[idx(i, j)] = s0 * (t0 * d0[idx(ix0, iy0)] + t1 * d0[idx(ix0, iy1)]) +
                s1 * (t0 * d0[idx(ix1, iy0)] + t1 * d0[idx(ix1, iy1)]);
        }
    }
    setBoundary(b, d);
}

fn project(u: []f32, v: []f32, p: []f32, div: []f32) void {
    const half_grid = 0.5 * @as(f32, @floatFromInt(N));

    var j: usize = 1;
    while (j <= N) : (j += 1) {
        var i: usize = 1;
        while (i <= N) : (i += 1) {
            div[idx(i, j)] = divergenceAt(u, v, i, j);
            p[idx(i, j)] = 0;
        }
    }

    setBoundary(0, div);
    setBoundary(0, p);
    linSolve(0, p, div, 1.0, 4.0);

    j = 1;
    while (j <= N) : (j += 1) {
        var i: usize = 1;
        while (i <= N) : (i += 1) {
            u[idx(i, j)] -= half_grid * (p[idx(i + 1, j)] - p[idx(i - 1, j)]);
            v[idx(i, j)] -= half_grid * (p[idx(i, j + 1)] - p[idx(i, j - 1)]);
        }
    }

    setBoundary(1, u);
    setBoundary(2, v);
}

fn stepFluid() void {
    diffuse(1, &vx_prev, &vx, VISCOSITY, DT);
    diffuse(2, &vy_prev, &vy, VISCOSITY, DT);
    project(vx_prev[0..], vy_prev[0..], &vx, &vy);

    advect(1, &vx, &vx_prev, &vx_prev, &vy_prev, DT);
    advect(2, &vy, &vy_prev, &vx_prev, &vy_prev, DT);
    project(vx[0..], vy[0..], &vx_prev, &vy_prev);

    diffuse(0, &density_prev, &density, DIFFUSION, DT);
    advect(0, &density, &density_prev, &vx, &vy, DT);

    total_density = 0;
    for (&density) |*d| {
        d.* *= FADE;
        if (d.* < MIN_DENSITY) d.* = 0;
        if (d.* > MAX_DENSITY) d.* = MAX_DENSITY;
        total_density += d.*;
    }
}

fn stepGame(dt: f32, input_x: f32, input_y: f32, reduced_motion: bool) void {
    updatePlayer(dt, input_x, input_y, reduced_motion);
    updateCamera();
    updateProximity();
    syncViews();
}

export fn gameUpdate(dt_raw: f32, input_x: f32, input_y: f32, reduced_motion_flag: u32) void {
    const dt = clamp(dt_raw, 0.0, 0.05);
    const reduced_motion = reduced_motion_flag != 0;
    stepGame(dt, input_x, input_y, reduced_motion);
    render_clock += dt;
    renderFrame();
}

export fn updateAndRender(dt_ms: f32) void {
    const dt = clamp(dt_ms / 1000.0, 0.0, 0.05);
    stepGame(dt, input_axis_x, input_axis_y, reducedMotionEnabled());
    render_clock += dt;
    renderFrame();
}

test "framebuffer dimensions and pointer contract" {
    try std.testing.expectEqual(@as(usize, Screen.width * Screen.height), FRAMEBUFFER_SIZE);
    try std.testing.expectEqual(@as(usize, Screen.width * Screen.height), getFramebufferLen());
    _ = getFramebufferPtr();
    _ = getPalettePtr();
}

test "rgb555 packs channels in gba bit layout" {
    try std.testing.expectEqual(@as(u16, 0x001f), rgb555(255, 0, 0));
    try std.testing.expectEqual(@as(u16, 0x03e0), rgb555(0, 255, 0));
    try std.testing.expectEqual(@as(u16, 0x7c00), rgb555(0, 0, 255));
}

test "set input clamps analog axes" {
    setInput(BUTTON_ACCEPT, 5.0, -6.0);
    try std.testing.expectEqual(@as(f32, 1.0), input_axis_x);
    try std.testing.expectEqual(@as(f32, -1.0), input_axis_y);
    try std.testing.expect((input_buttons & BUTTON_ACCEPT) != 0);
}

test "init clears game, effects, and fluid state" {
    init();
    addDensity(1, 1, 42.0);
    addVelocity(1, 1, 3.0, -2.0);
    gameUpdate(0.2, 1.0, 0.0, 0);

    init();

    try std.testing.expectEqual(@as(f32, 0), total_density);
    try std.testing.expectEqual(@as(f32, 0), density[idx(1, 1)]);
    try std.testing.expectEqual(@as(f32, 416), player.x);
    try std.testing.expect(getActivePortal() == -1);
}

test "density writes stay inside the active grid" {
    init();
    addDensity(1, 1, 10.0);
    addDensity(0, 1, 99.0);
    addDensity(N + 1, 1, 99.0);

    try std.testing.expectEqual(@as(f32, 10), density[idx(1, 1)]);
    try std.testing.expectEqual(@as(f32, 0), density[idx(0, 1)]);
    try std.testing.expectEqual(@as(f32, 0), density[idx(N + 1, 1)]);
}

test "projection divergence uses grid-normalized scale" {
    init();
    vx[idx(2, 1)] = 64.0;

    const expected = -32.0 / @as(f32, @floatFromInt(N));
    try std.testing.expectApproxEqAbs(expected, divergenceAt(vx[0..], vy[0..], 1, 1), 0.0001);
}

test "game update clamps movement inside world" {
    init();
    setPlayerPosition(1190.0, 790.0);
    gameUpdate(1.0, 1.0, 1.0, 1);

    try std.testing.expect(player.x <= World.width - 26.0);
    try std.testing.expect(player.y <= World.height - 26.0);
}

test "portal proximity activates and dismissal suppresses until exit" {
    init();
    setPlayerPosition(houses[0].x, houses[0].y);
    updateAndRender(16.0);

    try std.testing.expectEqual(@as(i32, 0), getActivePortal());
    dismissActivePortal();
    updateAndRender(16.0);
    try std.testing.expectEqual(@as(i32, -1), getActivePortal());

    setPlayerPosition(houses[0].x - TRIGGER_RADIUS - HYSTERESIS - 10.0, houses[0].y);
    updateAndRender(16.0);
    setPlayerPosition(houses[0].x, houses[0].y);
    updateAndRender(16.0);
    try std.testing.expectEqual(@as(i32, 0), getActivePortal());
}
