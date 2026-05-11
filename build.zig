const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{
        .default_target = .{
            .cpu_arch = .wasm32,
            .os_tag = .freestanding,
        },
    });

    const optimize = b.standardOptimizeOption(.{});

    const root_module = b.createModule(.{
        .root_source_file = b.path("fluid.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addExecutable(.{
        .name = "fluid",
        .root_module = root_module,
    });

    lib.entry = .disabled;
    lib.rdynamic = true;

    b.installArtifact(lib);
    b.installFile("index.html", "index.html");
    b.installFile("style.css", "style.css");
    b.installFile("game.js", "game.js");
    b.installFile("preloader.js", "preloader.js");

    const test_module = b.createModule(.{
        .root_source_file = b.path("fluid.zig"),
        .target = b.graph.host,
        .optimize = .Debug,
    });
    const tests = b.addTest(.{
        .root_module = test_module,
    });
    const run_tests = b.addRunArtifact(tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);
}
