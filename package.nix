{
  lib,
  stdenv,

  cmake,
  doxygen,

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
  version = "3.8";

  src = lib.cleanSource ./.;

  # With nix, CMAKE_INSTALL_* are absolute, so we must not prepend `${prefix}/`
  postPatch = ''
    substituteInPlace cmake/pkgconfig/simbody.pc.in --replace-fail "$""{prefix}/" ""
  '';

  nativeBuildInputs = [
    cmake
    doxygen
  ];

  buildInputs = [
    blas
    lapack
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
