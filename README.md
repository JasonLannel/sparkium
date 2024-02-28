# Sparks Renderer

## Advanced Computer Graphics: Final Project

### Description
This project is a neat path tracer, supporting
- basic features, including diffusive and specular materials;
- Bounding Volume Hierarchy with Surface Area Heuristic for
acceleration;
- Multiple importance sampling with power heuristic;
- Customed materials, including transmissive and Disney BSDF
materials;
- Customed textures, including color texture and normal map;
- Special visual effects, including motion blur, depth of field,
and alpha shadow;
- Anti-aliasing, with multi-sampling and FXAA;
- Simple implementation on homogeneous participant media;
- Classic scenes like Cornell Box and some customed scenes
for demo.

## Build The Project

Clone the repo (With submodules recursively):

```
git clone link-to-your-repo.git --recursive
```

or

```
git clone link-to-your-repo.git
cd sparkium
git submodule update --init --recursive
```

### Using Visual Studio 2022

Make sure that you installed Visual Studio with component `Desktop development with C++`.
If not, please open Visual Studio Installer and `Modify` your installation.

The repo contains configuration file `CMakeSettings.json`.
You only need to open the repo folder with Visual Studio on your local machine,
then Visual Studio should start CMake configuration process automatically.

When CMake configuration is successfully finished,
you should be able to select build type and debug target.
You should select `sparks.exe` as the debug target and run the program.

### Using Visual Studio Code

On Windows, you also should install Visual Studio with component `Desktop development with C++`.
We do this for the MSVC compiler, which has the best compatibility on Windows.
MinGW and MSYS2 are not recommended!!!

The repo also contained a vscode profile.
You can just open the repo folder on you local machine and configure CMake through settings on the bottom bar.
Your VSCode should have installed extensions `CMake Tools` and `C/C++` from Microsoft.

To run your program, you should select target `sparks` first. Then, you're recommended to use `Run->Start Debugging` and `Run->Run Without Debugging`, start from the menu bar.
