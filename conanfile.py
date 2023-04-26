from conan import ConanFile


class Recipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps", "VirtualRunEnv"

    def layout(self):
        self.folders.generators = "conan"

    def requirements(self):
        self.requires("fmt/9.1.0")
        self.requires("eigen/3.4.0")
        self.requires("shapelib/1.5.0")
        self.requires("boost/1.81.0")

    def build_requirements(self):
        self.test_requires("catch2/3.3.1")
