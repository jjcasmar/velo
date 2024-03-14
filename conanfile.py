import os

from conan.tools.cmake import CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.env import VirtualRunEnv

from conan import ConanFile

# Force Conan 2.0
required_conan_version = ">=2.0.0"


class Velo(ConanFile):
    name = "Velo"
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("entt/3.7.1")
        self.requires("spdlog/1.11.0")
        self.requires("tracy/0.10")
        self.requires("highfive/2.8.0")
        self.requires("pybind11/2.11.1")
    
    def configure(self):
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = True
        self.options["highfive"].with_xtensor = False
        self.options["highfive"].with_opencv = False

    def generate(self):
        cmake = CMakeDeps(self)
        cmake.generate()

        toolchain = CMakeToolchain(self)
        toolchain.user_presets_path = "ConanPresets.json"
        toolchain.generate()

        venv = VirtualRunEnv(self)
        venv.generate()

    def layout(self):
        cmake_layout(self)
