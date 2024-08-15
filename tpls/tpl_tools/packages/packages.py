class GenericPackage(object):

    def setDependencies(self, builder):
        return

    def configure_options(self, builder):
        return

    def configure(self, builder):
        return

    def build(self, builder):
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        return

    def install(self, builder):
        return

    def register(self, builder):
        raise NotImplementedError


class AutotoolsPackage(GenericPackage):

    def configure(self, builder):
        configure_options = ["./configure"]
        configure_options.extend(builder.get_options())
        configure_options.append("--prefix=" + builder.install_dir())
        builder.run_command(configure_options)

    def build(self, builder):
        builder.run_command(["make"], jobs_flag="-j", parallel=True)

    def install(self, builder):
        builder.run_command(["make", "install"])


class CMakePackage(GenericPackage):

    def configure(self, builder):
        cmake = builder._registry.get_executable("cmake")
        configure_options = [cmake, "-B", "build_tpl"]
        configure_options.extend(builder.get_options())
        configure_options.append("-DCMAKE_INSTALL_PREFIX=" + builder.install_dir())
        configure_options.append("-DCMAKE_BUILD_TYPE=Release")
        builder.run_command(configure_options)

    def build(self, builder):
        cmake = builder._registry.get_executable("cmake")
        builder.run_command(
            [cmake, "--build", "build_tpl"], jobs_flag="-j", parallel=True
        )

    def install(self, builder):
        cmake = builder._registry.get_executable("cmake")
        builder.run_command([cmake, "--install", "build_tpl"])
