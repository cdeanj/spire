from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext

setup(name = "c_hasher",
	ext_modules = cythonize("hashers/chash.pyx", include_dirs=['hashers']),
)