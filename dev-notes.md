# Notes

- ~~Set `BLA_VENDOR` to choose BLAS libraries when multiple BLAS's are installed~~
    - ~~Options at https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors~~

- CI updates/improvements
    - ~~Cancel concurrent runs~~
        ```
        concurrency:
          group: ${{ github.workflow }}-${{ github.ref }}
          cancel-in-progress: ${{ !contains(github.ref, 'simbody-')}}
        ```
    - Add single(?) CI run with oldest supported CMake
        ```
        - name: Install oldest supported CMake (test compatibility of CMakeLists.txt)
          uses: lukka/get-cmake@latest
          with:
            cmakeVersion: "~3.21"
        ```
    - Fix use of deprecated `set-output` https://github.blog/changelog/2022-10-11-github-actions-deprecating-save-state-and-set-output-commands/
    - Add ccache and cache-ing

- CMake visibility control
    - GCC/CLang, `-fvisibility=hidden` and `-fvisibility-inlines-hidden`
        - Add `SimTK_<LIBRARY>_EXPORT` define for GCC/Clang as `__attibute__((visibility("default")))`

- CMake RPATH setting?
    ```cmake
    if(NOT APPLE)
        set(CMAKE_INSTALL_RPATH $ORIGIN)
    endif()
    ```

- Ninja coloring: (From google/bloaty)
    ```
    # When using Ninja, compiler output won't be colorized without this.
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-fdiagnostics-color=always SUPPORTS_COLOR_ALWAYS)
    if(SUPPORTS_COLOR_ALWAYS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
    ```

- Component installs?
    ```cmake
    include(GNUInstallDirs)
    install(TARGETS Example
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            COMPONENT SomeProj_RunTime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            COMPONENT SomeProj_RunTime
            NAMELINK_COMPONENT SomeProj_Development
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            COMPONENT SomeProj_Development
        )
    ```

- Reproducible builds
    - https://blog.conan.io/2019/09/02/Deterministic-builds-with-C-C++.html
    - https://reproducible-builds.org/docs/
    - http://blog.llvm.org/2019/11/deterministic-builds-with-clang-and-lld.html
    - Diffoscope compare binary's
    - File ordering matters?
    - LTO
        - Creates random symbols?
        - Set SOURCE PROPERTY COMPILE_FLAGS `"-frandom-seed=0x${sha1sum of source file}"`
    - Macros
        - __DATE__ and __TIME__
            - MSVC: `add_link_options("/Brepro")`
            - GCC/Clang
                - Set env `SOURCE_DATE_EPOCH=<epoch_timestamp>`
                - `add_definitions(__DATE__="'<date>'")`, `add_definitions(__TIME__="'<time>'")`, and `add_compile_options(-Wno-builtin-macro-redefined)`
        - __FILE__


## Potential issues/PRs

- Add feature info's and summary?
    - https://cmake.org/cmake/help/latest/module/FeatureSummary.html
    - https://cmake.org/cmake/help/latest/module/FeatureSummary.html#command:add_feature_info
- Compiler flags
    - Current/past behavior overwrites manually specified flags
    - ~~Compiler flags specified in simbody aren't different from internal CMake platform/compiler defaults~~
        - This is true for GCC, Clang; not true for MSVC
    - Change debug info settings?
        - `-g1` for GCC, or `-gline-tables-only` for Clang in Release?
- Move adhoc tests behind build option? Not run by default, so no sense in building them by default
- VS Code debug configurations?
- Clang fmt, and analyzers (static, ASAN, UBSAN, etc)?
- Overlinking:
    - SimTKmath should only need to link to SimTKcommon? Math libs will be transitively linked due to public link of SimTKcommon. Same for Simbody
- `INSTALL_DOCS` isn't evenly/universally applied
    - Why doxygen separate/manually called instead of implied by `INSTALL_DOCS`?
- include-what-you-use analysis?
- Automatic release build artifact generation and built attestation?
    - Matrix:
        - RelWithDebInfo (or Release)
        - Shared and static
        - OS (Linux, MacOS, Windows)
        - Architecture?
            - ARM for every OS
            - x86 for every(? inc. MacOS?) OS
- Rewrite build instructions in README
- Make Simbody relocatable
    - simbody-visualizer hardcodes INSTALL_PREFIX

0. ~~NFC Cleanup~~
    - ~~Move `SIMBODY_(MAJOR|MINOR|PATCH)_VERSION` variables to project declaration and use `CMAKE_PROJECT_VERSION_(MAJOR|MINOR|PATCH)~~
        - ~~And same for subprojects and their versions~~
    - ~~Move `BUILD_USING_OTHER_LAPACK` up to top~~
    - ~~Remove text inside `else` and `endif` parentheses~~
    - ~~Mostly complete, just needs PR write up/summary~~
    - ~~Change set(VAR ${VAR} <items...>) style appends to correct list(APPEND VAR <items...>)~~
    - Add .git-blame-ignore-revs here?

1. ~~Add CMakePresets for configure, build, and test and update CI (Branch name: cmakepresets)~~
    - ~~Remove Travis and Appveyor config files. No longer in use afaik~~
        - ~~Is opensim actually using the AppVeyor nuget feed?~~
            - No; last mention of NuGet in opensim-core is from 2017/2018
        - ~~Grab any relevant config options from Travis file~~
            - ~~Coverage?~~
    - ~~Update GitHub actions CI to use presets~~
    - Helpful links:
        - https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html
        - https://martin-fieber.de/blog/cmake-presets/
    - ~~Update README?~~
    - ~~Confirm generality of adding BLA_VENDOR=Generic to build config?~~
    - Integrate with downstream builds? i.e. add presets for downstream?
    - ~~Remove parallel compilation flag from MSVC build flags (controlled in command line invocation)~~

2. Remove simultaneous static/shared builds and change namespace/versioning to modify the main/only set of libraries (i.e. not produce a copy) (Branch name: simplify-targets)
    - Add CMakePresets and documentation in README for how to replicate previous build configs (e.g. static and shared builds in single install, etc.)
        - Add "Static" configs? (Would enable dual static and shared installs when using multi-config buildsystems).
    - Announce changes on Simbody SimTK forum?
    - ~~Context/potentially cleaner methods to do so:~~
        - https://discourse.cmake.org/t/one-library-target-building-both-static-and-shared/3155
        - https://alexreinking.com/blog/building-a-dual-shared-and-static-library-with-cmake.html
        - https://stackoverflow.com/questions/2152077/is-it-possible-to-get-cmake-to-build-both-a-static-and-shared-library-at-the-sam
    - ~~Library versioning context~~
        - ~~`SOVERSION` target property may be  (but already in use?)~~
        - ~~`<$CONFIG>_POSTFIX` could replicate some functionality (the versioning suffix/postfix)~~
        - ~~static postfix using `ARCHIVE_OUTPUT_NAME`~~
    - ~~Update BUILD_TESTING docstring~~
    - ~~Add C/CXX_STANDARD(_REQUIRED)? properties for Simbody library targets~~
    - ~~`SimTK_SIMBODY_COPYRIGHT_YEARS` and `SimTK_SIMBODY_AUTHORS` compile definitions probably don't need to be project global?~~
    - ~~Move PlatformFiles target definition into the Platform CMakeLists.txt~~
        - ~~Blame to see why it's where it is and not in Platform~~
    - ~~Update SimbodyConfig.cmake.in~~ (and simbody.pc.in?)
    - ~~Add target headers via `target_sources` and `FILE_SET`s?~~
        - ~~Would need a version fence~~
    - ~~Rename <library>_INCLUDE_DIRS to something like <library>_BUILD_INCLUDE_DIRS to make~~
      ~~the BUILD_INTERFACE nature more explicit at time of use~~
    - ~~Every library needs to link(?) to depend on PlatformFiles?~~
        - ~~dirent.h not used? Confirmed.~~
            - ~~`PLATFORM_INCLUDE_DIRECTORIES` only set when `MSVC`, *but*, `#include <dirent.h>` only when ndef(_WIN32), and no other includes~~
    - ~~Test BUILD_SHARED_LIBS setting defaults for SIMBODY_BUILD_SHARED_LIBS~~

3. Finish changes needed to play nice as subproject 
    - Rename rest of `BUILD_` prefixed cache variables to use `SIMBODY_` prefix
    - BUILD_TESTING only if `Simbody_IS_TOP_LEVEL`
    - ~~BUILD_INTERFACE properties~~
    - ~~EXPORTS~~
    - May need to do in concert with 6
    - ~~Apply INSTALL_DOCS to other appropriate areas~~

4. Clean up build/install folder shenanigans
    - `CMAKE_INSTALL_DOCDIR`, ~~`BUILD_BINARY_DIR`, `EXECUTABLE_OUTPUT_PATH`, `LIBRARY_OUTPUT_PATH`,~~ `SIMBODY_VISUALIZER_ABS_INSTALL_DIR`
    - ~~`CMAKE_INSTALL_PREFIX` based modification of `CMAKE_INSTALL_LIBDIR`, `CMAKE_INSTALL_FULL_LIBDIR`. (Likely no longer necessary with modern CMake based on comments)~~

5. Add examples as tests (prevent bit/API-rot of examples)

6. Stop vendoring actual binary libraries
    - Make `WINDOWS_USE_EXTERNAL_LIBS` apply to glut as well?
    - Install BLAS and glut (freeglut?) on Windows
        - vcpkg has BLAS and freeglut packages; Windows Github Actions runners have vcpkg installed

- ~~Add configure time check to validate BLAS/LAPACK choice in BUILD_USING_OTHER_LAPACK~~ (branch name: trycompile-blas)
    - Successful check will print status messages during configure that read: 
      `-- Trying to compile with requested BLAS/LAPACK libraries - success.`
    - ```cmake
        message(CHECK_START "Trying to compile with requested BLAS/LAPACK libraries")
            
            # Try to link against the requested BLAS/LAPACK libraries 
            set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
            try_compile(OTHER_LAPACK_FUNCTIONAL
                SOURCE_FROM_CONTENT lapack_test.c
        "\
        #include \"SimTKlapack.h\"
        int main () {
            int* n = 2;
            int* stride = 1;
            float x[] = {1.0, 1.0};
            float y[] = {2.0, 3.0};

            float d = sdot_(n, x, stride, y, stride);

            return !(d == 5.0); // return zero if sdot_ worked
        }
        "
                CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${CMAKE_CURRENT_SOURCE_DIR}/SimTKcommon/include
                LINK_LIBRARIES ${LAPACK_BEING_USED})
            if(OTHER_LAPACK_FUNCTIONAL)
                message(CHECK_PASS "success.")
            else()
                message(CHECK_FAIL "failed.")
                message(FATAL_ERROR "The BLAS/LAPACK libraries requested using BUILD_USING_OTHER_LAPACK \
        did not compile successfully. Double-check that BUILD_USING_OTHER_LAPACK was specified correctly.")
            endif()
        ```
    - Add ` --debug-trycompile` option to debug failed checks that should work
        - Failed `try_compile`s will print the directory of the successful/failed try_compile project. `cd` to the directory and run `cmake --build .` to run the test build and see the build errors.
    - Testing:
        - These build configurations should succeed:
          `cmake -B build -S . --fresh` 
          `cmake -B build -S . --fresh -DBUILD_USING_OTHER_LAPACK="-lblas;-llapack"`
        - This build config should fail:
          `cmake -B build -S . --fresh -DBUILD_USING_OTHER_LAPACK="-lnotablas"`

## Downstream dependents/consumers to be aware of

- vcpkg => https://github.com/microsoft/vcpkg/blob/master/ports/simbody/portfile.cmake
    - Removes the need for patching
- Other recipes listed at https://repology.org/project/simbody/packages
- simbody conda-forge => https://github.com/conda-forge/simbody-feedstock
    - https://github.com/conda-forge/simbody-feedstock/blob/716dfa65bb306e308ae5408ea4184cc31347a3c6/recipe/meta.yaml#L78C9-L78C26 dirent not needed after PR#2
- opensim
- ~~Gazebo~~
    - Only Gazebo-classic depends on simbody, is maintained, but README directs users to the next generation of Gazebo (gazebo-sim?). Planned CMake changes would not affect use of simbody by gazebo-classic.
- MMB
    - Uses `findpackage`, would not be affected by planned CMake changes


# PRs
