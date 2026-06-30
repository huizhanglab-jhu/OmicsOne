from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

dummy_ext = Extension("omicsone_streamlit._binstub", sources=[])

class NoBuild(build_ext):
    def build_extensions(self):
        pass  # 不编译，只为标记 non-pure

setup(
    ext_modules=[dummy_ext],
    cmdclass={"build_ext": NoBuild},
)
