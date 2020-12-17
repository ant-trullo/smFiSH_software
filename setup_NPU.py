import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# ext_modules = [
#     Extension(
#         "hello",
#         ["hello.pyx"],
#         extra_compile_args=['-fopenmp'],
#         extra_link_args=['-fopenmp'],
#     )
# ]


setup(
    ext_modules=cythonize("NucsPileUpUtility.pyx", annotate=True),
    include_dirs=[numpy.get_include()]
)

# python3 setup.py build_ext --inplace
