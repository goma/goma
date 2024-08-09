import os


class Config(object):
    def __init__(self):
        self.environment = {}

    def set_environment_variable(self, variable, value):
        self.environment[variable] = value

    def append_environment_variable(self, variable, value):
        if variable in self.environment.keys():
            if type(self.environment[variable]) != list:
                oldval = self.environment[variable]
                self.environment[variable] = [oldval]
            self.environment[variable].append(value)
        else:
            self.environment[variable] = value

    def write_config(self, file, shell="bash"):
        with open(file, "w") as f:
            if shell == "bash":
                for k in self.environment.keys():
                    if type(self.environment[k]) == list:
                        f.write("export {}=".format(k))
                        for item in self.environment[k]:
                            f.write(item)
                            f.write(":")
                        f.write("${}".format(k))
                        f.write("\n")
                    else:
                        f.write("export {}=".format(k))
                        f.write("{}".format(self.environment[k]))
                        f.write("\n")
            elif shell == "fish":
                for k in self.environment.keys():
                    if type(self.environment[k]) == list:
                        f.write("set -x {} ".format(k))
                        for item in self.environment[k]:
                            f.write(item)
                            f.write(" ")
                        f.write("${}".format(k))
                        f.write("\n")
                    else:
                        f.write("export {}=".format(k))
                        f.write("{}".format(self.environment[k]))
                        f.write("\n")


class Registry(object):
    def __init__(self):
        self.executables = {}
        self.package_locations = {}
        self.environment = os.environ.copy()
        self.config = Config()

    def register_executable(self, exe_path):
        if os.path.isfile(exe_path):
            name = os.path.basename(exe_path)
            self.executables[name] = exe_path
            print("Registering {}:{}".format(name, exe_path))
        else:
            raise ValueError("exe not a file {}".format(exe_path))

    def register_package(self, package_name, package_path):
        if os.path.exists(package_path):
            self.package_locations[package_name] = package_path

    def get_package(self, package_name):
        return self.package_locations[package_name]

    def get_executable(self, exe):
        return self.executables[exe]

    def get_environment(self):
        return self.environment

    def set_environment_variable(self, variable, value):
        self.environment[variable] = value
        self.config.set_environment_variable(variable, value)

    def append_environment_variable(self, variable, value):
        if variable in self.environment:
            self.environment[variable] += os.pathsep + value
        else:
            self.environment[variable] = value
        self.config.append_environment_variable(variable, value)
