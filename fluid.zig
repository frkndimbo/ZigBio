const std = @import("std");

const Screen = struct {
    const width: f32 = 240.0;
    const height: f32 = 160.0;
};

const World = struct {
    const width: f32 = 900.0;
    const height: f32 = 620.0;
};

const N: usize = 64;
const SIZE: usize = (N + 2) * (N + 2);
const ITER: usize = 4;
const HOUSE_COUNT: usize = 5;
const RIPPLE_COUNT: usize = 72;
const PARTICLE_COUNT: usize = 72;

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
    .{ .x = 442.0, .y = 330.0, .radius = TRIGGER_RADIUS },
    .{ .x = 590.0, .y = 300.0, .radius = TRIGGER_RADIUS },
    .{ .x = 720.0, .y = 395.0, .radius = TRIGGER_RADIUS },
    .{ .x = 380.0, .y = 472.0, .radius = TRIGGER_RADIUS },
    .{ .x = 650.0, .y = 535.0, .radius = TRIGGER_RADIUS },
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

fn clamp(value: f32, min: f32, max: f32) f32 {
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

export fn init() void {
    player = .{
        .x = 520.0,
        .y = 360.0,
        .vx = 0.0,
        .vy = 0.0,
        .dir = 0.0,
        .moving = 0.0,
        .step_clock = 0.0,
        .ripple_clock = 0.0,
    };
    camera = .{ 400.0, 280.0 };
    active_house = -1;
    dismissed_house = -1;
    resetFluid();
    resetEffects();
    syncHouseView();
    syncViews();
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

export fn getActiveHouse() i32 {
    return active_house;
}

export fn dismissActiveHouse() void {
    dismissed_house = active_house;
    active_house = -1;
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
    const was_moving = player.moving > 0.5;
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
    player.ripple_clock += dt;

    if (moving and player.ripple_clock > 0.11) {
        player.ripple_clock = 0.0;
        spawnRipple(player.x, player.y, input_len, reduced_motion);
        disturbFluid(player.x, player.y, player.vx, player.vy);
    }

    if (was_moving and !moving) {
        spawnSplash(player.x, player.y + 8.0, reduced_motion);
    }
}

fn updateCamera() void {
    const target_x = clamp(player.x - Screen.width / 2.0, 0.0, World.width - Screen.width);
    const target_y = clamp(player.y - Screen.height / 2.0, 0.0, World.height - Screen.height);
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

export fn gameUpdate(dt_raw: f32, input_x: f32, input_y: f32, reduced_motion_flag: u32) void {
    const dt = clamp(dt_raw, 0.0, 0.05);
    const reduced_motion = reduced_motion_flag != 0;
    updatePlayer(dt, input_x, input_y, reduced_motion);
    updateCamera();
    updateEffects(dt);
    stepFluid();
    updateProximity();
    syncViews();
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

test "init clears game, effects, and fluid state" {
    init();
    addDensity(1, 1, 42.0);
    addVelocity(1, 1, 3.0, -2.0);
    gameUpdate(0.2, 1.0, 0.0, 0);

    init();

    try std.testing.expectEqual(@as(f32, 0), total_density);
    try std.testing.expectEqual(@as(f32, 0), density[idx(1, 1)]);
    try std.testing.expectEqual(@as(f32, 520), player.x);
    try std.testing.expect(active_house == -1);
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

    try std.testing.expectApproxEqAbs(@as(f32, -0.5), divergenceAt(vx[0..], vy[0..], 1, 1), 0.0001);
}

test "game update clamps movement inside world" {
    init();
    setPlayerPosition(1190.0, 790.0);
    gameUpdate(1.0, 1.0, 1.0, 1);

    try std.testing.expect(player.x <= World.width - 26.0);
    try std.testing.expect(player.y <= World.height - 26.0);
}

test "house proximity activates and dismissal suppresses until exit" {
    init();
    setPlayerPosition(houses[0].x, houses[0].y);
    gameUpdate(0.016, 0.0, 0.0, 1);

    try std.testing.expectEqual(@as(i32, 0), active_house);
    dismissActiveHouse();
    gameUpdate(0.016, 0.0, 0.0, 1);
    try std.testing.expectEqual(@as(i32, -1), active_house);

    setPlayerPosition(houses[0].x - TRIGGER_RADIUS - HYSTERESIS - 10.0, houses[0].y);
    gameUpdate(0.016, 0.0, 0.0, 1);
    setPlayerPosition(houses[0].x, houses[0].y);
    gameUpdate(0.016, 0.0, 0.0, 1);
    try std.testing.expectEqual(@as(i32, 0), active_house);
}
