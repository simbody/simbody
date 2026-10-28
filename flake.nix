{
  description = "High-performance C++ multibody dynamics/physics library for simulating articulated biomechanical and mechanical systems like vehicles, robots, and the human skeleton";

  inputs = {
    nixpkgs.url = "https://channels.nixos.org/nixpkgs-unstable/nixexprs.tar.zst";
  };

  outputs = inputs: {
    packages = builtins.mapAttrs (system: pkgs: {
      simbody = pkgs.callPackage ./package.nix { };

      default = inputs.self.packages.${system}.simbody;
    }) inputs.nixpkgs.legacyPackages;
  };
}
