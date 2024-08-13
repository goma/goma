import os
import subprocess
import tarfile
import sys
from tpl_tools import utils


def mkdir_p(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def print_and_log(line, file):
    sys.stdout.write(line)
    file.write(line)


def print_and_log_err(line, file):
    sys.stderr.write(line)
    file.write(line)


class Builder(object):

    def __init__(
        self,
        package,
        jobs,
        download_dir,
        extract_dir,
        install_dir,
        logger,
        registry,
        build_shared,
    ):
        self._dependencies = []
        self._configure_options = []
        self._package = package
        self.build_shared = build_shared
        self._download_dir = download_dir
        self._extract_dir = extract_dir
        self._extracted_folder = None
        self._install_dir = install_dir
        self.logger = logger
        self._jobs = jobs
        self._registry = registry
        return

    def set_dependency(self, dep):
        self._dependencies.append(dep)

    def get_dependencies(self):
        return list(set(self._dependencies))

    def add_option(self, opt_string):
        self._configure_options.append(opt_string)

    def get_options(self):
        return self._configure_options

    def download(self):
        mkdir_p(self._download_dir)
        url = self._package.url
        sha256 = self._package.sha256
        filename = os.path.join(self._download_dir, self._package.filename)
        self.logger.log("Downloading file: {}".format(filename))
        success = utils.download_file(url, filename, sha256)
        if success:
            self.logger.log("Successfully downloaded: {}".format(filename))
        return success

    def extract(self):
        mkdir_p(self._extract_dir)
        filename = os.path.join(self._download_dir, self._package.filename)
        self.logger.log("Extracting {}".format(self._package.name))
        self.logger.log("Extracting from {}".format(filename))
        tf = tarfile.open(filename)
        tf.extractall(os.path.join(self._extract_dir))
        self._extracted_folder = tf.getnames()[0]
        self.logger.log(
            "Extracted {} to {}".format(
                self._package.name,
                os.path.join(self._extract_dir, self._extracted_folder),
            )
        )

    def build(self):
        self._package.set_environment(self)
        self._package.configure_options(self)
        if hasattr(self._package, "configure"):
            self.logger.log("Configuring {}".format(self._package.name))
            self._package.configure(self)
        if hasattr(self._package, "build"):
            self.logger.log("Building {}".format(self._package.name))
            self._package.build(self)

    def install(self):
        if hasattr(self._package, "install"):
            self.logger.log("Installing {}".format(self._package.name))
            self._package.install(self)

    def register(self):
        if hasattr(self._package, "register"):
            self.logger.log("Registering:", self._package.name)
            self._package.register(self)

    def install_dir(self):
        return os.path.join(
            self._install_dir, self._package.name + "-" + self._package.version
        )

    def run_command(
        self, command_list, jobs_flag=None, parallel=False, env=None, command_id=None
    ):
        command = command_list
        if parallel and jobs_flag:
            if jobs_flag[-1] == "=":
                jobs_flag = jobs_flag + str(self._jobs)
                command.append(jobs_flag)
            else:
                command.append(jobs_flag)
                command.append(str(self._jobs))
        self.logger.log("Running command: {}".format(" ".join(command)))
        cwd = os.path.join(self._extract_dir, self._extracted_folder)
        self.logger.log("In directory: {}".format(cwd))
        log_directory = os.path.join(self._install_dir, "logs")
        mkdir_p(log_directory)
        logfile = self._package.name + "-" + self._package.version
        if command_id:
            logfile += str(command_id)
        logfile += ".log"
        logfile = os.path.join(log_directory, logfile)
        self.logger.log("log available at : {}".format(logfile))
        with open(logfile, "a") as lf:
            runenv = self.env
            if env:
                runenv = env
            cmd = subprocess.Popen(command, cwd=cwd, env=runenv, stderr=lf, stdout=lf)
        cmd.wait()
        if cmd.returncode != 0:
            with open(logfile, "r") as f:
                for line in f:
                    sys.stderr.write(line)
            raise RuntimeError("Command resulted in a non-zero error code")

    def check(self, skip_named_install=False):
        self.logger.log("Checking install of {}".format(self._package.name))
        install_dir = self.install_dir()
        if skip_named_install:
            install_dir = self._install_dir
        if hasattr(self._package, "executables"):
            bindir = os.path.join(install_dir, "bin")
            for e in self._package.executables:
                exe = os.path.join(bindir, e)
                self.logger.log("Checking for {}".format(exe))
                if not os.path.isfile(exe):
                    self.logger.log("File not found: ", exe)
                    return False
                self.logger.log("Found executable: {}".format(exe))
        if hasattr(self._package, "libraries"):
            for lib in self._package.libraries:
                found = False
                prefix = ["", "lib"]
                directories = ["", "lib", "lib64"]
                extension = utils.get_library_extension(self.build_shared)
                self.logger.log("Checking for {}".format(lib))
                for d in directories:
                    for p in prefix:
                        checkpath = install_dir
                        if d != "":
                            checkpath = os.path.join(checkpath, d)
                        if os.path.exists(os.path.join(checkpath, p + lib + extension)):
                            found = True
                if not found:
                    self.logger.log(
                        "Failed to find library: {}".format("lib" + lib + extension)
                    )
                    return False
                self.logger.log("Found library: {}".format(lib))
        if hasattr(self._package, "includes"):
            for header in self._package.includes:
                directories = ["", "include"]
                found = False
                for d in directories:
                    checkpath = install_dir
                    if d != "":
                        checkpath = os.path.join(checkpath, d)
                    if os.path.exists(os.path.join(checkpath, header)):
                        found = True
                if not found:
                    self.logger.log("Failed to find include: {}".format(header))
                    return False
                self.logger.log("Found include: {}".format(header))
        return True
