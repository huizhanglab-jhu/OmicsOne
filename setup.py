from setuptools import setup
from setuptools_rust import Binding, RustExtension


setup(
    rust_extensions=[
        RustExtension(
            "omicsone.utils._spearmanr",
            "packages/rust_spearmanr/Cargo.toml",
            binding=Binding.PyO3,
            debug=False,
        )
    ],
    zip_safe=False,
)
