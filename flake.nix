{
  description = "Flake for R development";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    systems.url = "github:nix-systems/default";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
    nix-gl-host = {
      url = "github:numtide/nix-gl-host";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { nixpkgs, flake-utils, nix-gl-host, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };

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
            cudaPackages.cuda_nvcc
            cudaPackages.cuda_cudart
            cudaPackages.libcurand
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
            # run your cuda application with 'nixglhost ./your_app'
            nix-gl-host.defaultPackage.x86_64-linux
          ];
        };
      });
}
