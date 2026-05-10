const std = @import("std");

// Grid configuration - using 64 for better performance and visibility
const N: usize = 64;
const SIZE: usize = (N + 2) * (N + 2);
const ITER: usize = 4;

// Fluid properties - tuned for visible, long-lasting fluid
const DIFFUSION: f32 = 0.001;
const VISCOSITY: f32 = 0.0001;
const DT: f32 = 0.2;
const FADE: f32 = 0.992; // Slower fade for longer trails
const MIN_DENSITY: f32 = 0.01;
const MAX_DENSITY: f32 = 1000.0;

// Memory buffers
var density: [SIZE]f32 = undefined;
var density_prev: [SIZE]f32 = undefined;
var vx: [SIZE]f32 = undefined;
var vy: [SIZE]f32 = undefined;
var vx_prev: [SIZE]f32 = undefined;
var vy_prev: [SIZE]f32 = undefined;

// Exposed status metric for UI and tests.
var total_density: f32 = 0;

export fn getTotalDensity() f32 {
    return total_density;
}

// Helper function to get 1D index from 2D coordinates
fn idx(x: usize, y: usize) usize {
    return x + (N + 2) * y;
}

// Initialize all buffers to zero
export fn init() void {
    @memset(&density, 0);
    @memset(&density_prev, 0);
    @memset(&vx, 0);
    @memset(&vy, 0);
    @memset(&vx_prev, 0);
    @memset(&vy_prev, 0);
    total_density = 0;
}

// Get pointer to density array for JavaScript
export fn getDensityPtr() [*]u8 {
    return @ptrCast(&density);
}

// Get grid size
export fn getGridSize() usize {
    return N;
}

// Add density at specific grid cell
export fn addDensity(x: usize, y: usize, amount: f32) void {
    if (x >= 1 and x <= N and y >= 1 and y <= N) {
        density[idx(x, y)] += amount;
    }
}

// Add velocity at specific grid cell
export fn addVelocity(x: usize, y: usize, amountX: f32, amountY: f32) void {
    if (x >= 1 and x <= N and y >= 1 and y <= N) {
        vx[idx(x, y)] += amountX;
        vy[idx(x, y)] += amountY;
    }
}

// Set boundary conditions
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

// Linear solver using Gauss-Seidel relaxation
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

// Diffusion step
fn diffuse(b: i32, x: []f32, x0: []f32, diff: f32, dt: f32) void {
    const a = dt * diff * @as(f32, @floatFromInt(N * N));
    linSolve(b, x, x0, a, 1.0 + 4.0 * a);
}

fn divergenceAt(u: []const f32, v: []const f32, i: usize, j: usize) f32 {
    const inv_grid = 1.0 / @as(f32, @floatFromInt(N));
    return -0.5 * inv_grid * (u[idx(i + 1, j)] - u[idx(i - 1, j)] + v[idx(i, j + 1)] - v[idx(i, j - 1)]);
}

// Advection step (backtrace method)
fn advect(b: i32, d: []f32, d0: []f32, u: []f32, v: []f32, dt: f32) void {
    const dt0 = dt * @as(f32, @floatFromInt(N));

    var j: usize = 1;
    while (j <= N) : (j += 1) {
        var i: usize = 1;
        while (i <= N) : (i += 1) {
            var x = @as(f32, @floatFromInt(i)) - dt0 * u[idx(i, j)];
            var y = @as(f32, @floatFromInt(j)) - dt0 * v[idx(i, j)];

            if (x < 0.5) x = 0.5;
            if (x > @as(f32, @floatFromInt(N)) + 0.5) x = @as(f32, @floatFromInt(N)) + 0.5;

            if (y < 0.5) y = 0.5;
            if (y > @as(f32, @floatFromInt(N)) + 0.5) y = @as(f32, @floatFromInt(N)) + 0.5;

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

// Projection step for mass conservation (incompressibility)
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

// Main simulation step
export fn step() void {
    // Velocity step
    diffuse(1, &vx_prev, &vx, VISCOSITY, DT);
    diffuse(2, &vy_prev, &vy, VISCOSITY, DT);

    project(vx_prev[0..], vy_prev[0..], &vx, &vy);

    advect(1, &vx, &vx_prev, &vx_prev, &vy_prev, DT);
    advect(2, &vy, &vy_prev, &vx_prev, &vy_prev, DT);

    project(vx[0..], vy[0..], &vx_prev, &vy_prev);

    // Density step
    diffuse(0, &density_prev, &density, DIFFUSION, DT);
    advect(0, &density, &density_prev, &vx, &vy, DT);

    // Apply fade and track total
    total_density = 0;
    for (&density) |*d| {
        d.* *= FADE;
        if (d.* < MIN_DENSITY) d.* = 0;
        if (d.* > MAX_DENSITY) d.* = MAX_DENSITY;
        total_density += d.*;
    }
}

test "init clears simulation state" {
    init();
    addDensity(1, 1, 42.0);
    addVelocity(1, 1, 3.0, -2.0);

    init();

    try std.testing.expectEqual(@as(f32, 0), total_density);
    try std.testing.expectEqual(@as(f32, 0), density[idx(1, 1)]);
    try std.testing.expectEqual(@as(f32, 0), vx[idx(1, 1)]);
    try std.testing.expectEqual(@as(f32, 0), vy[idx(1, 1)]);
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

test "step keeps density finite and clamped" {
    init();
    addDensity(N / 2, N / 2, MAX_DENSITY * 4.0);
    addVelocity(N / 2, N / 2, 8.0, -6.0);

    step();

    try std.testing.expect(total_density > 0);
    for (density) |value| {
        try std.testing.expect(!std.math.isNan(value));
        try std.testing.expect(value <= MAX_DENSITY);
    }
}
