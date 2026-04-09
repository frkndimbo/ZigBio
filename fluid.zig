const std = @import("std");

// Grid configuration
const N: usize = 128;
const SIZE: usize = (N + 2) * (N + 2);
const ITER: usize = 4;

// Fluid properties
const DIFFUSION: f32 = 0.00001;
const VISCOSITY: f32 = 0.000001;
const DT: f32 = 0.1;

// Memory buffers
var density: [SIZE]f32 = undefined;
var density_prev: [SIZE]f32 = undefined;
var vx: [SIZE]f32 = undefined;
var vy: [SIZE]f32 = undefined;
var vx_prev: [SIZE]f32 = undefined;
var vy_prev: [SIZE]f32 = undefined;

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

            const i0: usize = @intFromFloat(x);
            const i1 = i0 + 1;
            const j0: usize = @intFromFloat(y);
            const j1 = j0 + 1;

            const s1 = x - @as(f32, @floatFromInt(i0));
            const s0 = 1.0 - s1;
            const t1 = y - @as(f32, @floatFromInt(j0));
            const t0 = 1.0 - t1;

            d[idx(i, j)] = s0 * (t0 * d0[idx(i0, j0)] + t1 * d0[idx(i0, j1)]) +
                           s1 * (t0 * d0[idx(i1, j0)] + t1 * d0[idx(i1, j1)]);
        }
    }
    setBoundary(b, d);
}

// Projection step for mass conservation (incompressibility)
fn project(u: []f32, v: []f32, p: []f32, div: []f32) void {
    var j: usize = 1;
    while (j <= N) : (j += 1) {
        var i: usize = 1;
        while (i <= N) : (i += 1) {
            div[idx(i, j)] = -0.5 * (@as(f32, @floatFromInt(N))) * (u[idx(i + 1, j)] - u[idx(i - 1, j)] + v[idx(i, j + 1)] - v[idx(i, j - 1)]);
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
            u[idx(i, j)] -= 0.5 * (@as(f32, @floatFromInt(N))) * (p[idx(i + 1, j)] - p[idx(i - 1, j)]);
            v[idx(i, j)] -= 0.5 * (@as(f32, @floatFromInt(N))) * (p[idx(i, j + 1)] - p[idx(i, j - 1)]);
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
}
