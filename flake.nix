{
  description = "Flake for R development";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.systems.url = "github:nix-systems/default";
  inputs.flake-utils = {
    url = "github:numtide/flake-utils";
    inputs.systems.follows = "systems";
  };

  outputs = { nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        microsimulation = pkgs.rPackages.buildRPackage {
          name = "microsimulation";
          src = pkgs.fetchFromGitHub {
            owner = "mclements";
            repo = "microsimulation";
            rev = "27211dd9941ddda6e8a494bbb6462c281650fdd4";
            hash = "sha256-pYo5kRq1kQNDbrqWfwU9Ydk2kKSRNhm5sUfrycVhWds=";
          };
          propagatedBuildInputs = with pkgs.rPackages; [ ascii survival RcppArmadillo ];
        };
        prostata = pkgs.rPackages.buildRPackage {
          name = "prostata";
          src = ./.;
          propagatedBuildInputs = [ microsimulation ];
        };

        vscDebugger = pkgs.rPackages.buildRPackage {
          name = "vscDebugger";
          src = pkgs.fetchFromGitHub {
            owner = "ManuelHentschel";
            repo = "vscDebugger";
            rev = "v0.5.8";
            hash = "sha256-+eq+V33eb9r+7193YTOzJFzmFqZxRm8UjAq+8C95qZ0=";
          };
          propagatedBuildInputs = with pkgs.rPackages; [ jsonlite R6 ];
        };

      in {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [ 
            gtest
            zstd
            xz
            bzip2
            zlib
            icu
            rPackages.ascii
            rPackages.survival
            rPackages.RcppArmadillo
            microsimulation
            R
            rPackages.devtools
            rPackages.dplyr
            rPackages.minqa
            vscDebugger
            prostata #In R, run: devtools::load_all('prostata')"
          ];
        };
      });
}
