{
  lib,
  stdenv,

  cmake,

  blas,
  lapack,
  freeglut,
  libglvnd,
  libGLU,
  libxmu,
  libxi,
}:

stdenv.mkDerivation (finalAttrs: {
  pname = "simbody";
  version = "3.9.0";

  src = lib.cleanSource ./.;

  nativeBuildInputs = [
    cmake
  ];

  propagatedBuildInputs = [
    blas
    lapack
  ];

  buildInputs = [
    freeglut
    libglvnd
    libGLU
    libxmu
    libxi
  ];

  cmakeFlags = [
    (lib.cmakeBool "BUILD_SHARED_LIBS" (!stdenv.hostPlatform.isStatic))
    (lib.cmakeBool "BUILD_TESTING" finalAttrs.finalPackage.doCheck)
  ];

  doCheck = true;

  meta = {
    description = "High-performance C++ multibody dynamics/physics library for simulating articulated biomechanical and mechanical systems like vehicles, robots, and the human skeleton";
    homepage = "https://github.com/simbody/simbody";
    license = lib.licenses.asl20;
    platforms = lib.platforms.unix ++ lib.platforms.windows;
  };
})
