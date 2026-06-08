{
  description = "ERP para fábricas de costura";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-parts.url = "github:hercules-ci/flake-parts";
  };

  outputs = inputs @ {flake-parts, ...}:
    flake-parts.lib.mkFlake {inherit inputs;} {
      systems = ["x86_64-linux"];

      perSystem = {
        config,
        pkgs,
        lib,
        system,
        ...
      }: {
        devShells.default = pkgs.mkShell {
          buildInputs = [pkgs.ffmpeg_8];
          PKG_CONFIG_PATH = "${pkgs.ffmpeg_8}/lib/pkgconfig";
        };
      };
    };
}
