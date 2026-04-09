const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });

    const optimize = b.standardOptimizeOption(.{});

    const lib = b.addExecutable(.{
        .name = "fluid",
        .root_source_file = b.path("fluid.zig"),
        .target = target,
        .optimize = optimize,
    });

    lib.entry = .disabled;
    lib.rdynamic = true;

    b.installArtifact(lib);
}
