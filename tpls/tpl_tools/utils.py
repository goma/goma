import urllib
import urllib.request
import shutil
import os
import hashlib
import sys
import subprocess
import re


# https://stackoverflow.com/questions/22058048/hashing-a-file-in-python
def sha256sum(filename, buffer_size=65536):
    sha256 = hashlib.sha256()
    with open(filename, "rb", buffering=0) as f:
        while True:
            data = f.read(buffer_size)
            if not data:
                break
            sha256.update(data)
    return sha256.hexdigest()


def check_gcc_clang_version(cc):
    version = subprocess.check_output([cc, "--version"]).decode("utf-8")
    vn = re.search(r"(\d+\.\d+\.\d+)", version).group(1)
    if "gcc" in version:
        return "gcc", tuple(map(int, vn.split(".")))
    if "clang" in version:
        return "clang", tuple(map(int, vn.split(".")))
    return None, vn


VERBOSE = True


def topological_sort(dep_graph):
    # Create a copy of dependencies and initialize in-degrees
    in_degree = {pkg: 0 for pkg in dep_graph}
    for deps in dep_graph.values():
        for dep in deps:
            in_degree[dep] += 1

    # Start with packages that have no dependencies
    order = []
    queue = [pkg for pkg, deg in in_degree.items() if deg == 0]

    while queue:
        pkg = queue.pop(0)
        order.append(pkg)
        for neighbor in dep_graph.get(pkg, []):
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    if len(order) != len(dep_graph):
        raise Exception("Cyclic dependency detected!")

    return order


def set_verbosity(verbosity):
    global VERBOSE
    VERBOSE = verbosity


def print_if_verbose(*args):
    if VERBOSE:
        print(*args)


def get_library_extension(shared):
    if not shared:
        return ".a"
    if sys.platform == "darwin":
        return ".dylib"
    return ".so"


def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
    return None


def download_file(url, filename, sha256=None, verify=True):
    skip_download = False
    if os.path.exists(filename):
        print("")
        if os.path.isfile(filename):
            print_if_verbose("File already exists: {}".format(filename))
            skip_download = True

    if not skip_download:
        if verify:
            with urllib.request.urlopen(url) as req, open(filename, "wb") as f:
                shutil.copyfileobj(req, f)
        else:
            import ssl

            context = ssl.create_default_context()
            context.check_hostname = False
            context.verify_mode = ssl.CERT_NONE
            with urllib.request.urlopen(url, context=context) as req:
                with open(filename, "wb") as f:
                    shutil.copyfileobj(req, f)

    if os.path.isfile(filename):
        if not sha256 is None:
            sum = sha256sum(filename)
            if sha256 != sum:
                print_if_verbose("SHA mismatch {}".format(filename))
                print_if_verbose("Expected: {}".format(sha256))
                print_if_verbose("Actual: {}".format(sum))
                raise ValueError(sha256 + " != " + sum)
        return True
    else:
        return False


class PrintAndFileLogger(object):
    def __init__(self, filepath):
        self.file = open(filepath, "a")

    def log(self, *args):
        print(*args)
        print(*args, file=self.file)

    def __del__(self):
        self.file.close()


class PrintLogger(object):
    def __init__(self):
        return

    def log(self, *args):
        print(*args)


def check_for_x11(extract_dir, cc):
    comp = "cc"
    if cc:
        comp = cc
    with open(os.path.join(extract_dir, "x11test.c"), "w") as f:
        f.write("int main() { return 0; }\n")
    command = [cc, "-lX11", "x11test.c"]
    result = subprocess.Popen(command, cwd=extract_dir)
    result.wait()
    os.remove(os.path.join(extract_dir, "x11test.c"))
    return result.returncode == 0
